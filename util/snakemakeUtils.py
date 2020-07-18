import json
import time
import os
import errno
import csv
import re

def dict_merge(dct, merge_dct):
    """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], dict)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]


# writes the config dict into a file in the output_dir with time stamp, need to be called after the config is populated:
def write_config_file(config_dict, out_dir_key="out_dir"):
    # create output dir if it is not exist:
    date_str = time.strftime("%b-%d-%Y-%H-%M", time.localtime())
    if os.path.basename(os.getcwd())==config_dict[out_dir_key]:
        config_file_name = "config_{}.json".format(date_str)
    else:   
        os.makedirs(config_dict[out_dir_key],exist_ok=True)
        config_file_name = os.path.join(config_dict[out_dir_key], "config_%s.json" % date_str)
        
    with open(config_file_name, 'w') as pipeline_config_file_handle:
        json.dump(config_dict, pipeline_config_file_handle, indent=4)

def job_name(job_name):
    """Job name in LSF should be truncated since is limited in size"""
    return job_name[:15] if len(job_name) > 15 else job_name


def is_key(config, key):
    return key in config and config[key].strip()

class StringMatcher:
    """
        StringMatcher: defines a string matching class in a list which is exact and case insensitive.
        returns a list with all the matches or the [] for no match
        Usage:
            header_matcher = StringMatcher(['sample', 'lane', 'file'])
            header_matcher['SaMplE']
            >>['sample']
            header_matcher['Sam']
            >>[]
    """
    def __init__(self,strlist):
        self.strlist = strlist
    def __getitem__(self, key):
        # returns all (case-insensitive) matches
        re_str=r"(?i)\b{}\b".format(key)
        return [re.match(re_str, item).group() for item in self.strlist if re.match(re_str, item) is not None]
    def __str__(self):
        print("StringMatcher for the list (1st 10 items):")
        return '{}'.format(self.strlist[:10])
    def __repr__(self):
        return self.__str__()
    

class SampleInfoReader(object):
    """
        Sample config structure is hierarchical
            sample_config['sample']['col_i']
            Examples:
                sample_config['sample']['read1'] # fastq, read1 file 
                sample_config['sample']['fastq'] # fastq file 
                sample_config['sample']['bam'] # bam file ...
    """
    # TBD: continue developing the methods for hierarchical sample config dict
    
    @staticmethod
    def old_sample_info_reader(file_path, delimiter=',', header=False, paired=False, file_type='fastq'):
        """
        The order of cols in the file needs to be: sample, fastq1[, fastq2]
        :param delimiter: the delimiter (',' or '\t')
        :param header: if there is a header or not
        :param paired: if it is paired or not
        :return: None
        """
        samples_config = {}
        num_cols = 3 if paired else 2
        with open(file_path, 'r', newline='') as samples_info_handle:
            csv_reader = csv.reader(samples_info_handle, delimiter=delimiter)
            if header:
                next(csv_reader)
            for (i, row) in enumerate(csv_reader):
                row = [elem.strip() for elem in row]
                if len(row) != num_cols:
                    raise Exception(
                        "There are %s columns in line %s of the samples_info_file instead of %s(name, path1[,path2]), "
                        "don't forget to use %s as a delimiter!" % (str(len(row)), str(i + 1), str(num_cols), delimiter))
                if not os.path.isfile(row[1]) or (paired and not os.path.isfile(row[2])):
                    raise Exception(
                        "one of the paths in line %s is not a file, don't forget to remove the header!" % str(i + 1))
                if paired:
                    samples_config[row[0]] = {file_type: [row[1], row[2]]}
                else:
                    samples_config[row[0]] = {file_type: row[1]}
        return samples_config

    @staticmethod
    def sample_path_reader(filename, delimiter='\t'):
        """
            reads a table from a file with a header
            the output is exported to a dict with sample as a key and path as a value, 
            therefore sample,path must be in the columns.
            the format of the ouptut is
            sample_config['sample']=path            
        """
        with open(filename, newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=delimiter)
            header_matcher = StringMatcher(reader.fieldnames)
            # find the sample row
            sample_match = header_matcher['Sample']
            if len(sample_match)!=1:
                raise Exception('Error: One of the columns in the input file must include the word <sample>!')
            sample_str=sample_match[0]
            # find the path row
            path_match = header_matcher['path']
            if len(sample_match)!=1:
                raise Exception('Error: One of the columns in the input file must include the word <path>!')
            path_str=path_match[0]
            # create the dict
            sample_dict={}
            for row in reader:
                #print('{}: {}'.format(row[sample_str], row[path_str]))
                sample_dict[row[sample_str]]=row[path_str]
        return sample_dict
        
    @staticmethod
    def sample_table_reader(filename, delimiter='\t', key_name='sample', col_names=['path'], opt_col_names=[]):
        """
            Reads a table from a file with a header
            One of the columns will be used as the main key and must match the argument key_name.
            The output is exported to a dict with key_name as the main key. 
            The other given columns (required or optional) will be used in a nested dictionary as follows:
            D[key1][col1]=val11
            D[key1][col2]=val12, etc.
        """
        with open(filename, newline='') as csvfile:
            reader = csv.DictReader((row for row in csvfile if not row.startswith('#')), delimiter=delimiter)
            #lower_case_fieldnames = [item.lower() for item in reader.fieldnames]
            header_matcher = StringMatcher(reader.fieldnames)
            # find the column for the key (usually sample)
            key_match = header_matcher[key_name.lower()]
            if len(key_match)!=1:
                raise Exception('Error: One of the columns in the input file must include the word: %s'%format(key_name))        
            key_idx=key_match[0]
            # find the obligatory columns
            col_strs=[]
            for col in col_names:
                # find the matched col name
                col_match = header_matcher[col.lower()]
                if len(col_match)!=1:
                    raise Exception('Error: One of the columns in the input file must include the word %s!'%format(col))
                col_strs.append(col_match[0])
            # add optional columns (if any)
            for col in opt_col_names:
                # find the matched col name
                col_match = header_matcher[col.lower()]
                if len(col_match)==1:
                    col_strs.append(col_match[0])
            # create the dict
            fileinfo_dict={}
            for row in reader:
                fileinfo_dict[row[key_idx]]={col.lower():row[col].strip() for col in col_strs}
            return fileinfo_dict

