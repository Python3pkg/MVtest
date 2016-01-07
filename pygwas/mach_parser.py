import pygwas
from pygwas.data_parser import DataParser
from parsed_locus import ParsedLocus
from . import sys_call
from . import ExitIf
import sys
from exceptions import TooManyAlleles
from exceptions import TooFewAlleles
import gzip
import numpy
from exceptions import InvalidSelection
import os

from pheno_covar import PhenoCovar

__copyright__ = "Eric Torstenson"
__license__ = "GPL3.0"
#     This file is part of pyGWAS.
#
#     pyGWAS is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     pyGWAS is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

# Fake enum
class Encoding(object):
    #: Currently there is only one way to interpret these values
    Dosage = 0

encoding = Encoding.Dosage

class Parser(DataParser):
    """Parse IMPUTE style output.

    Due to the nature of the mach data format, we must load the data first
    into member before we can begin analyzing it. Due to the massive amount
    of data, SNPs are loaded in in chunks.

    ISSUES:
        * Currently, we will not be filtering on individuals except by explicit removal
        * We are assuming that each gzip archive contains all data associated with the loci contained within (i.e. there won't
            be separate files with different subjects inside) (( Todd email jan-9-2015))
        * There is no reason to process regions in any order. I'm thinking we'll have a master file and then indices into that
            file and task count to facilitate "parallel" execution

        * There is no place to store RSID from the output that I've seen (Minimac output generated by Ben Zhang)


    """

    #: Extension for the info file
    info_ext = "info.gz"
    #: Extension for the dosage file
    dosage_ext = "dose.gz"
    #: Number of loci to parse at a time (larger stride requires more memory)
    chunk_stride = 50000
    #: rsquared threshold for analysis (obtained from the mach output itself)
    min_rsquared = 0.3

    def __init__(self, archive_list, info_files=[]):
        """Initialize the structure with the family details file and the list of archives to be parsed

        :param archive_list:    list of gzipped files to be considered
        :info_files:            (optional) list of corresponding info files. If the extensions are correct, we can
                                derive the info filename
        :return:

        This function assumes that all dosage files have the same sample order (as with the output from minimac)
        """
        smallest = -1
        self.family_details = None      #  we'll use the smallest file to determine the sample order

        infos = []                      # temporary to build up info_files if they aren't provided
        idx = 0

        if not DataParser.compressed_pedigree:
            if Parser.dosage_ext[-3:] == ".gz":
                Parser.dosage_ext = Parser.dosage_ext[0:-3]
            if Parser.info_ext[-3:] == ".gz":
                Parser.info_ext = Parser.info_ext[0:-3]

        for file in archive_list:
            try:
                s = os.stat(file).st_size

                if self.family_details is None or smallest > s:
                    self.family_details = file
                    smallest = s

                if len(info_files) == 0:
                    info_file = file.replace(Parser.dosage_ext, Parser.info_ext)
                    pygwas.ExitIf("Info file not found, %s" % (info_file), not os.path.exists(info_file))
                    pygwas.ExitIf("Info and sample files appear to be same. Is the gen_ext invalid? (%s)" % info_file, info_file == file)
                    infos.append(info_file)
                idx+=1
            except:
                pygwas.ExitIf("Archive file not found, %s" % (file), True)

        self.archives = archive_list        # This is only the list of files to be processed
        if len(info_files) > 0:
            self.info_files = info_files
        else:
            self.info_files = infos
        self.current_file = self.archives[0]            # This will be used to record the opened file used for parsing
        self.current_info = self.info_files[0]            # This will be used to record the info file associated with quality of SNPs

        self.chunk = 0
        self.file_index = 0

        assert len(self.info_files) == len(self.archives)


    def ReportConfiguration(self, file):
        """Report the configuration details for logging purposes.

        :param file: Destination for report details
        :return: None
        """
        global encodingpar
        print >> file, pygwas.BuildReportLine("MACH_ARCHIVES", "%s")
        idx = 0
        for arch in self.archives[0:]:
            print >> file, pygwas.BuildReportLine("", "%s:%s" % (self.archives[idx], self.info_files[idx]))
            idx += 1
        print >> file, pygwas.BuildReportLine("ENCODING", ["Dosage", "Genotype"][encoding])


    def load_family_details(self, pheno_covar):
        """Load contents from the .fam file, updating the pheno_covar with \
            family ids found.

        :param pheno_covar: Phenotype/covariate object
        :return: None
        """
        self.file_index = 0
        mask_components = []        # 1s indicate an individual is to be masked out

        #print >> sys.stderr, "Loading Family Details"

        file = self.family_details
        if DataParser.compressed_pedigree:
            data, serr = sys_call('gunzip -c %s | wc -l' % (file))
            self.line_count = int(data[0].strip().split(" ")[0])
            iddata, serr = sys_call('gunzip -c %s | cut -f 1' % (file))
        else:
            data, serr = sys_call('wc -l %s' % (file))
            self.line_count = int(data[0].strip().split(" ")[0])
            iddata, serr = sys_call('cat %s | cut -f 1' % (file))

        ids_observed = set()
        for line in iddata:
            indid = line.strip().split()[0]
            indid = ":".join(indid.split("->"))

            ExitIf("Duplicate ID found in dose file: %s" % (indid), indid in ids_observed)
            ids_observed.add(indid)

            if DataParser.valid_indid(indid):
                mask_components.append(0)
                pheno_covar.add_subject(indid, PhenoCovar.missing_encoding, PhenoCovar.missing_encoding)
            else:
                mask_components.append(1)

        self.ind_mask = numpy.array(mask_components) == 1
        self.ind_count = self.ind_mask.shape[0]
        pheno_covar.freeze_subjects()

        #print >> sys.stderr, "Family Details Loaded"

    def openfile(self, filename):
        if DataParser.compressed_pedigree:
            return gzip.open(filename, 'rb')
        return open(filename, 'r')

    def parse_genotypes(self, lb, ub):
        """Extracts a fraction of the file (current chunk of loci) loading
        the genotypes into memoery.

        :param lb: Lower bound of the current chunk
        :param ub: Upper bound of the current chunk
        :return: Dosage dosages for current chunk

        """
        file = self.openfile(self.current_file)
        words = file.readline().strip().split()[lb:ub]
        word_count = len(words)
        idx =0

        if word_count > 0:
            dosages = numpy.empty((self.ind_count, word_count), dtype='|S5')
            while word_count > 1:
                dosages[idx] = numpy.array(words)
                idx += 1
                line = file.readline()
                words = line.strip().split()[lb:ub]
                word_count = len(words)
        else:
            raise EOFError

        return dosages

    def load_genotypes(self):
        """Actually loads the first chunk of genotype data into memory due to \
        the individual oriented format of MACH data.

        Due to the fragmented approach to data loading necessary to avoid
        running out of RAM, this function will initialize the data structures
        with the first chunk of loci and prepare it for otherwise normal
        iteration.

        Also, because the parser can be assigned more than one .gen file to
        read from, it will automatically move to the next file when the
        first is exhausted.

        """

        lb = self.chunk * Parser.chunk_stride + 2
        ub = (self.chunk + 1) * Parser.chunk_stride + 2

        buff = None

        self.current_file = self.archives[self.file_index]
        self.info_file = self.info_files[self.file_index]

        while buff is None:
            try:
                buff = self.parse_genotypes(lb, ub)
            except EOFError:
                buff = None
                if self.file_index < (len(self.archives) - 1):
                    self.file_index += 1
                    self.chunk = 0
                    lb = self.chunk * Parser.chunk_stride + 2
                    ub = (self.chunk + 1) * Parser.chunk_stride + 2
                    self.current_file = self.archives[self.file_index]
                    self.info_file = self.info_files[self.file_index]
                else:
                    raise StopIteration

        # Numpy's usecols don't prevent it from loading entire file, which is too big considering ours are 60+ gigs
        #buff = numpy.loadtxt(self.current_file, usecols=range(lb, ub), dtype=str)
        self.dosages = numpy.transpose(buff)

        file = self.openfile(self.info_file)
        file.readline()     # drop header

        lindex = 0
        while lindex < lb - 2:
            file.readline()
            lindex += 1

        self.markers = []
        self.rsids = []
        self.locus_count= 0
        self.maf = []
        self.alleles = []
        self.rsquared = []

        while lindex < (ub - 2):
            words = file.readline().strip().split()
            if len(words) > 0:
                loc, al2, al1, freq1, maf, avgcall,rsq = words[0:7]
                marker = loc.split(":")[0:2]
                marker[0]=int(marker[0])
                self.markers.append(marker)
                self.maf.append(float(maf))
                self.alleles.append([al1, al2])
                self.rsquared.append(float(rsq))
                lindex += 1
            else:
                break

        if self.dosages.shape[0] != len(self.markers):
            print >> sys.stderr, "What is going on? I have ", self.dosages.shape[0], "dosages per individual and ", len(self.markers), self.markers

        self.chunk += 1
        self.marker_count = len(self.markers)

        #print >> sys.stderr, "Genotypes loaded"

    def get_effa_freq(self, genotypes):
        """Returns the frequency of the effect allele"""
        return numpy.mean(numpy.array(genotypes)/2)


    def populate_iteration(self, iteration):
        """Parse genotypes from the file and iteration with relevant marker \
            details.

        :param iteration: ParseLocus object which is returned per iteration
        :return: True indicates current locus is valid.

        StopIteration is thrown if the marker reaches the end of the file or
        the valid genomic region for analysis.

        This function will force a load of the next chunk when necessary.
        """
        global encoding

        # We never have to worry about iteration exceptions...that is taken care of in load_genotypes
        cur_idx = iteration.cur_idx
        if cur_idx < 0 or cur_idx >= self.marker_count:
            self.load_genotypes()
            iteration.cur_idx = 0
            cur_idx = 0

        iteration.chr = self.markers[cur_idx][0]
        iteration.pos = int(self.markers[cur_idx][1])

        if DataParser.boundary.TestBoundary(iteration.chr, iteration.pos, iteration.rsid) and self.rsquared[cur_idx] >= Parser.min_rsquared:
            iteration.major_allele, iteration.minor_allele = self.alleles[cur_idx]
            iteration.genotype_data = numpy.ma.masked_array(self.dosages[cur_idx].astype(numpy.float), self.ind_mask).compressed()
            iteration._maf = numpy.mean(iteration.genotype_data/2)
            iteration.allele_count2 = (iteration.genotype_data.shape[0] * 4.0 - numpy.sum(iteration.genotype_data))

            return iteration.maf >= DataParser.min_maf and iteration.maf <= DataParser.max_maf
        return False


    def __iter__(self):
        """Reset the file and begin iteration"""

        return ParsedLocus(self)