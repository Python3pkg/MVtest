#!/usr/bin/env python
import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../")
    sys.path.insert(0, "../../") # For mvtest
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import numpy
import os
import mvtest
from libgwas.boundary import BoundaryCheck
from libgwas import mach_parser
from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
import libgwas.standardizer
import gzip
import unittest
numpy.random.seed(1337)

base_freq = [0.95, 0.75, 0.7,0.8, 0.65, 0.7, 0.85, 0.7, 0.7, 0.3]


class TestBase(unittest.TestCase):
    def setUp(self):
        self.allele_1 = list("AAACCCGGGTCGTGTATACC")
        self.allele_2 = list("CGTGTATACCAAACCCGGGT")

        self.WriteTestFiles()

        self.phenotypes = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        self.sex = [1,2,1,1,2,2,1,1,2,2,2,1]
        self.chrpos_encoding = mach_parser.Parser.chrpos_encoding
        self.dosage_ext = mach_parser.Parser.dosage_ext
        self.info_ext = mach_parser.Parser.info_ext
        self.chrom = BoundaryCheck.chrom
        self.boundary = DataParser.boundary
        DataParser.boundary = BoundaryCheck()
        self.min_maf = DataParser.min_maf
        self.max_maf = DataParser.max_maf
        self.snp_miss_tol = DataParser.snp_miss_tol
        self.ind_miss_tol = DataParser.ind_miss_tol
        DataParser.ind_exclusions = []
        self.sex_as_covar = PhenoCovar.sex_as_covariate
        self.has_sex = DataParser.has_sex
        self.has_pheno = DataParser.has_pheno
        self.has_parents = DataParser.has_parents
        self.has_fid = DataParser.has_fid
        self.has_liability = DataParser.has_liability
        self.encoding = mach_parser.encoding
        self.compression = DataParser.compressed_pedigree
        DataParser.compressed_pedigree = True
        #self.chunk_stride = mach_parser.chunk_stride
        self.standardizer = libgwas.standardizer.get_standardizer()
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)

    def tearDown(self):
        os.remove(self.gen_file)
        os.remove(self.gen_file2)
        os.remove(self.uncmp_1)
        os.remove(self.uncmp_2)
        os.remove(self.info_file1)
        os.remove(self.info_file2)
        os.remove(self.info_ucmp1)
        os.remove(self.info_ucmp2)
        os.remove(self.mach_file)
        os.remove(self.pheno_covar)

        mach_parser.Parser.dosage_ext = self.dosage_ext
        mach_parser.Parser.info_ext = self.info_ext
        mach_parser.Parser.chrpos_encoding = self.chrpos_encoding

        BoundaryCheck.chrom  = self.chrom
        DataParser.boundary  = self.boundary
        DataParser.min_maf   = self.min_maf
        DataParser.max_maf   = self.max_maf
        DataParser.snp_miss_tol  = self.snp_miss_tol
        DataParser.ind_miss_tol  = self.ind_miss_tol
        DataParser.ind_exclusions = []
        DataParser.has_sex    = self.has_sex
        DataParser.has_pheno  = self.has_pheno
        DataParser.has_fid    = self.has_fid
        DataParser.has_liability = self.has_liability
        DataParser.has_parents = self.has_parents
        PhenoCovar.sex_as_covariate = self.sex_as_covar
        DataParser.compressed_pedigree = self.compression
        #mach_parser.chunk_stride = self.chunk_stride
        libgwas.standardizer.set_standardizer(self.standardizer)


    def WriteTestFiles(self, prefix = "__test_imputed"):

        self.ind_ids = ["ID0001->FAM001",
                        "ID0002->FAM002",
                        "ID0003->FAM003",
                        "ID0004->FAM004",
                        "ID0005->FAM005",
                        "ID0006->FAM006",
                        "ID0007->FAM007",
                        "ID0008->FAM008",
                        "ID0009->FAM009",
                        "ID0010->FAM010",
                        "ID0011->FAM011",
                        "ID0012->FAM012"]

        self.gen_file = "%s.dose.gz" % (prefix)
        self.gen_file2 = "%s-2.dose.gz" % (prefix)
        self.info_file1 = "%s.info.gz" % (prefix)
        self.info_file2 = "%s-2.info.gz" % (prefix)
        self.uncmp_1 = "%s.dose" % (prefix)
        self.uncmp_2 = "%s-2.dose" % (prefix)
        self.info_ucmp1 = "%s.info" % (prefix)
        self.info_ucmp2 = "%s-2.info" % (prefix)
        self.mach_file = "%s.mach" % (prefix)

        self.pheno_covar = "__phenocovar.txt"
        with open(self.pheno_covar, "w") as file:
            file.write("""FID\tIID\tSEX\tAGE\tBMI
ID0001\tFAM001\t1\t30\t28.54
ID0002\tFAM002\t1\t33\t30.10
ID0003\tFAM003\t2\t28\t24.00
ID0004\tFAM004\t2\t40\t29.21
ID0005\tFAM005\t1\t50\t31.23
ID0006\tFAM006\t1\t30\t29.54
ID0007\tFAM007\t1\t33\t33.10
ID0008\tFAM008\t2\t28\t27.00
ID0009\tFAM009\t2\t40\t27.21
ID0010\tFAM010\t1\t50\t30.23
ID0011\tFAM011\t2\t28\t24.00
ID0012\tFAM012\t2\t40\t29.21
""")

        with open(self.mach_file, "w") as file:
            print("%s.dose.gz %s.info.gz" % (prefix, prefix), file=file)
            print("%s-2.dose.gz %s-2.info.gz" % (prefix, prefix), file=file)

        gen_file = gzip.open(self.gen_file, 'wb')
        uncmp_file = open(self.uncmp_1, 'w')
        idx = 0
        self.dosage_encoding = numpy.zeros((20, 12))
        self.positions = []
        self.mafs = numpy.zeros(len(base_freq) * 2)

        info_file = open(self.info_file1, 'w')
        print("snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0", file=info_file)

        self.chroms = [ int(x) for x in ['1'] * 7 + ['2'] * 7 + ['3'] * 6]
        self.positions = [1012, 1020, 1026, 1032, 1100, 1137, 1149] * 2 + [1012, 1020, 1026, 1032, 1100, 1137]
        self.alleles = [list(numpy.random.choice(['A','C','G','T'], 2, replace=False)) for x in range(0, 20)]
        idx = 0

        mafs = numpy.zeros((10))
        dosages = numpy.zeros((12, 10))
        for ind in self.ind_ids:
            f = numpy.random.normal(base_freq, scale=0.1)
            f[f>=1.0] = 0.99
            maf = 1.0 - f
            AA = f * f
            Aa = 2 * f * maf
            aa = maf * maf
            dosages[idx] = Aa + 2*AA
            mafs += dosages[idx] / 2
            print("\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            ), file=gen_file)
            print("\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            ), file=uncmp_file)
            idx += 1
        self.mafs[0:10] = mafs/10
        self.dosage_encoding[0:10,:] = numpy.transpose(dosages)
        gen_file.close()
        uncmp_file.close()
        info_file = gzip.open(self.info_file1, 'wb')
        info_ufile = open(self.info_ucmp1, 'w')
        print("SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2", file=info_file)
        print("SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2", file=info_ufile)
        for idx in range(0, 10):
            print("\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ]), file=info_file)
            print("\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ]), file=info_ufile)
        info_file.close()
        info_ufile.close()


        gen_file = gzip.open(self.gen_file2, 'wb')
        uncmp_file = open(self.uncmp_2, 'w')

        idx = 0
        mafs = numpy.zeros((10))
        dosages = numpy.zeros((12, 10))
        for ind in self.ind_ids:
            f = numpy.random.normal(base_freq, scale=0.1)
            f[f>=1.0] = 0.99
            maf = 1.0 - f
            mafs += maf
            AA = f * f
            Aa = 2 * f * maf
            aa = maf * maf
            dosages[idx] = Aa + 2*aa
            print("\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            ), file=gen_file)
            print("\t".join([
                ind,
                "DOSE"] +
                ["%.3f" % x for x in dosages[idx]]
            ), file=uncmp_file)
            idx += 1
        self.mafs[10:] = mafs/10
        self.dosage_encoding[10:,:] = numpy.transpose(dosages)

        gen_file.close()

        info_file = gzip.open(self.info_file2, 'wb')
        info_cfile = open(self.info_ucmp2, 'w')
        print("SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2", file=info_file)
        print("SNP\tAl1\tAl2\tFreq1\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tdose2", file=info_cfile)
        for idx in range(10, 20):
            print("\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ]), file=info_file)
            print("\t".join([
                "%s:%d" % (self.chroms[idx],self.positions[idx]),
                self.allele_1[idx],
                self.allele_2[idx],
                str(1.0-self.mafs[idx]),
                str(self.mafs[idx]),
                '0.99912',
                '0.8',
                "\t".join(['-'] * 6)
            ]), file=info_cfile)
        info_cfile.close()
        info_file.close()

        self.mach_parser = mach_parser.Parser([self.gen_file])


class TestMachCmdLine(TestBase):
    def testMachCmdLineBasics(self):
        cmds = "--mach %s --pheno %s --pheno-names BMI" % \
               (self.mach_file, self.pheno_covar)
        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        self.assertEqual("PhenoCovar", vars.__class__.__name__)
        self.assertEqual(2, len(dataset.archives))
        self.assertEqual([self.gen_file, self.gen_file2], dataset.archives)
        self.assertEqual([self.info_file1, self.info_file2], dataset.info_files)

        self.assertEqual(12, dataset.ind_count)

    def testMachCmdLineChrPos(self):
        cmds = "--mach %s --pheno %s --pheno-names BMI --mach-chrpos --chr 1" % \
               (self.mach_file, self.pheno_covar)
        app = mvtest.MVTestApplication()

        # No exception because we can use chromosomes with chrpos
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        cmds ="--mach %s --pheno %s --pheno-names BMI --chr 1 --from-bp 1 --to-bp 10000 --mach-chrpos" % \
               (self.mach_file, self.pheno_covar)
        app = mvtest.MVTestApplication()

        # No exception because we can use chromosomes with chrpos
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        # You can't use positional boundaries without giving it a chromosome
        cmds ="--mach %s --pheno %s --pheno-names BMI --from-bp 1 --to-bp 10000 --mach-chrpos" % \
               (self.mach_file, self.pheno_covar)
        app = mvtest.MVTestApplication()
        with self.assertRaises(SystemExit) as cm:
            dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual(cm.exception.code, 1)

    def testMachCmdLineNoChrPos(self):
        cmds = "--mach %s --pheno %s --pheno-names BMI --chr 1" % \
               (self.mach_file, self.pheno_covar)
        app = mvtest.MVTestApplication()
        # This is illegal, since we can't use chromosomes along with default
        # MACH
        with self.assertRaises(SystemExit) as cm:
            dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual(cm.exception.code, 1)

        cmds = "--mach %s --pheno %s --pheno-names BMI --chr 1 --from-bp 1 --to-bp 10000" % \
               (self.mach_file, self.pheno_covar)
        app = mvtest.MVTestApplication()
        # This is illegal, since we can't use chromosomes along with default
        # MACH
        with self.assertRaises(SystemExit) as cm:
            dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual(cm.exception.code, 1)

    def testMachCmdLineMultiplePhenotypes(self):
        cmds = "--mach %s --pheno %s --pheno-names AGE,BMI" % \
               (self.mach_file, self.pheno_covar)
        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        index = 0
        for var in vars:
            index += 1
        self.assertEqual(2, index)



if __name__ == "__main__":
    unittest.main()
