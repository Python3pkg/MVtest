#!/usr/bin/env python

import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../../")
    sys.path.insert(0, "../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import unittest
import numpy
from . import test_transped_parser
from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
from libgwas.transposed_pedigree_parser import Parser as TransposedPedigreeParser
import libgwas.standardizer
from meanvar.mvstandardizer import Standardizer
import meanvar.mvstandardizer



class TestTPedStandardization(test_transped_parser.TestBase):
    def setUp(self):
        super(TestTPedStandardization, self).setUp()
        self.standardizer = libgwas.standardizer.get_standardizer()
        libgwas.standardizer.set_standardizer(Standardizer)



    def tearDown(self):
        super(TestTPedStandardization, self).tearDown()
        libgwas.standardizer.set_standardizer(self.standardizer)

    def test_tped_standardization(self):
        DataParser.has_sex = True
        DataParser.has_pheno = True
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        nonmissing = numpy.empty(pc.phenotype_data[0].shape, dtype=numpy.bool)
        nonmissing[:] = True
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)

        raw_pheno = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        raw_cov   = [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1]

        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))

            for i in range(0, len(raw_pheno)):
                self.assertAlmostEqual(raw_pheno[i], y[i])
                self.assertAlmostEqual(raw_cov[i], c[0][i])

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        pc.do_standardize_variables = True
        libgwas.standardizer.set_standardizer(Standardizer)

        std_pheno = [-1.61601695, -0.73455316,  1.02837442, -0.4407319 , 0.73455316, 1.02837442,
                     -1.61601695, -0.73455316,  1.02837442, -0.4407319 , 0.73455316, 1.02837442]
        std_cov   = [-0.70710678, -0.70710678,  1.41421356,  1.41421356, -0.70710678, -0.70710678,
                     -0.70710678, -0.70710678,  1.41421356,  1.41421356, -0.70710678, -0.70710678]
        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))

            for i in range(0, len(std_pheno)):
                self.assertAlmostEqual(std_pheno[i], y[i])
                self.assertAlmostEqual(std_cov[i], c[0][i])

    def test_tped_standardization2(self):
        DataParser.has_sex = True
        DataParser.has_pheno = True
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        nonmissing = numpy.empty(pc.phenotype_data[0].shape, dtype=numpy.bool)
        nonmissing[:] = True
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)

        raw_pheno = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        raw_cov   = [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1]

        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))

            for i in range(0, len(raw_pheno)):
                self.assertAlmostEqual(raw_pheno[i], y[i])
                self.assertAlmostEqual(raw_cov[i], c[0][i])

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        pc.do_standardize_variables = True
        libgwas.standardizer.set_standardizer(Standardizer)

        std_pheno = [-1.61601695, -0.73455316,  1.02837442, -0.4407319 , 0.73455316, 1.02837442,
                     -1.61601695, -0.73455316,  1.02837442, -0.4407319 , 0.73455316, 1.02837442]
        std_cov   = [-0.70710678, -0.70710678,  1.41421356,  1.41421356, -0.70710678, -0.70710678,
                     -0.70710678, -0.70710678,  1.41421356,  1.41421356, -0.70710678, -0.70710678]
        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))

            for i in range(0, len(std_pheno)):
                self.assertAlmostEqual(std_pheno[i], y[i])
                self.assertAlmostEqual(std_cov[i], c[0][i])

    def test_tped_standardization_w_missing2(self):
        PhenoCovar.sex_as_covariate = True

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        nonmissing = numpy.empty(pc.phenotype_data[0].shape, dtype=numpy.bool)
        nonmissing[:] = True
        nonmissing[0] = False
        nonmissing[1] = False
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)

        raw_pheno = [1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5, 0.9, 1.0]
        raw_cov   = [2, 2, 1, 1, 1, 1, 2, 2, 1, 1]

        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))
            for i in range(0, len(raw_pheno)):
                self.assertAlmostEqual(raw_pheno[i], y[i])
                self.assertAlmostEqual(raw_cov[i], c[0][i])

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        pc.do_standardize_variables = True
        libgwas.standardizer.set_standardizer(Standardizer)

        std_pheno = [1.02837442, -0.4407319 , 0.73455316, 1.02837442,
                     -1.61601695, -0.73455316,  1.02837442, -0.4407319 , 0.73455316, 1.02837442]
        std_cov   = [1.41421356,  1.41421356, -0.70710678, -0.70710678,
                     -0.70710678, -0.70710678,  1.41421356,  1.41421356, -0.70710678, -0.70710678]
        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))

            for i in range(0, len(std_pheno)):
                self.assertAlmostEqual(std_pheno[i], y[i])
                self.assertAlmostEqual(std_cov[i], c[0][i])


    def test_tped_standardization_w_missing1(self):
        PhenoCovar.sex_as_covariate = True
        DataParser.ind_exclusions = ["11:11", "12:12"]

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        nonmissing = numpy.empty(pc.phenotype_data[0].shape, dtype=numpy.bool)
        nonmissing[:] = True
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)

        raw_pheno = [0.1, 0.4, 1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5]
        raw_cov   = [1, 1, 2, 2, 1, 1, 1, 1, 2, 2]

        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))
            for i in range(0, len(raw_pheno)):
                self.assertAlmostEqual(raw_pheno[i], y[i])
                self.assertAlmostEqual(raw_cov[i], c[0][i])

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        pc.do_standardize_variables = True
        libgwas.standardizer.set_standardizer(Standardizer)

        std_pheno = [-1.43314068, -0.55570761,  1.19915853, -0.26322992,  0.90668084,
                        1.19915853, -1.43314068, -0.55570761,  1.19915853, -0.26322992]
        std_cov   = [-0.81649658, -0.81649658,  1.22474487,  1.22474487, -0.81649658,
                        -0.81649658, -0.81649658, -0.81649658,  1.22474487,  1.22474487]
        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))

            for i in range(0, len(std_pheno)):
                self.assertAlmostEqual(std_pheno[i], y[i])
                self.assertAlmostEqual(std_cov[i], c[0][i])

    def test_tped_standardization_w_dbl_missing(self):
        PhenoCovar.sex_as_covariate = True
        DataParser.ind_exclusions = ["11:11", "12:12"]

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()
        nonmissing = numpy.empty(pc.phenotype_data[0].shape, dtype=numpy.bool)
        nonmissing[:] = True
        nonmissing[0] = False
        nonmissing[1] = False
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)

        raw_pheno = [1.0, 0.5, 0.9, 1.0, 0.1, 0.4, 1.0, 0.5]
        raw_cov   = [2, 2, 1, 1, 1, 1, 2, 2]


        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))
            for i in range(0, len(raw_pheno)):
                self.assertAlmostEqual(raw_pheno[i], y[i])
                self.assertAlmostEqual(raw_cov[i], c[0][i])

        pc = PhenoCovar()
        ped_parser = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        ped_parser.load_tfam(pc)
        ped_parser.load_genotypes()

        pc.do_standardize_variables = True
        libgwas.standardizer.set_standardizer(Standardizer)

        std_pheno = [ 1.19915853, -0.26322992,  0.90668084,
                        1.19915853, -1.43314068, -0.55570761,  1.19915853, -0.26322992]
        std_cov   = [ 1.22474487,  1.22474487, -0.81649658,
                        -0.81649658, -0.81649658, -0.81649658,  1.22474487,  1.22474487]
        test_var  = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        for pheno in pc:
            (y, c, total_nonmissing) = pheno.get_variables(numpy.invert(nonmissing))


            for i in range(0, len(std_pheno)):
                self.assertAlmostEqual(std_pheno[i], y[i])
                self.assertAlmostEqual(std_cov[i], c[0][i])


if __name__ == "__main__":
    unittest.main()
