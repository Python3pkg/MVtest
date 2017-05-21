import sys
import math
__author__ = 'Eric Torstenson'
__version__ = "0.9.1"

__copyright__ = "Copyright (C) 2015 Todd Edwards, Chun Li and Eric Torstenson"
__license__ = "GPL3.0"
#     This file is part of MVtest.
#
#     MVtest is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     MVtest is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with MVtest.  If not, see <http://www.gnu.org/licenses/>.


class MVResult(object):
    """Result associated with a single locus/phenotype execution

    """

    def __init__(self, chr, pos, rsid, maj, min, eff_alcount, non_miss_count, p_mvtest, ph_label, beta_values , pvalues, stderrors, maf, covar_labels=[], lm=-1,runtime=-1):
        #: Chromosome
        self.chr = chr
        #: BP position
        self.pos = pos
        #: RSID
        self.rsid = rsid
        #: Major allele (A,C,G,T, etc)
        self.maj_allele = maj
        #: Minor allele
        self.min_allele = min
        #: Total count of effect alleles
        self.eff_alcount = eff_alcount
        #: non missing count
        self.non_miss = non_miss_count
        #: mvtest's pvalue
        self.p_mvtest = p_mvtest
        #: list of beta values
        self.betas = beta_values
        #: list of beta pvalues
        self.beta_pvalues = pvalues
        #: list of std errors
        self.beta_stderr = stderrors
        #: current phenotype label
        self.ph_label   = ph_label
        #: minor allele frequency
        self.maf = maf
        #: number of seconds analysis took to complete
        self.runtime = runtime
        #: LM
        self.lmpv = lm
        #: Covariate labels used for analysis
        self.covar_labels = covar_labels


    @property
    def p_variance(self):
        return self.beta_pvalues[3]

    def print_header(self, f=sys.stdout, verbose=False):
        """Prints header to f (will write header based on verbose)

        :param f: stream to print output
        :param verbose: print all data or only the most important parts?
        """
        self.var_count = 2 + len(self.covar_labels)

        if verbose:
            header = [
                "Chr",
                "Pos",
                "RSID",
                "Phenotype",
                "N",
                "Ref_allele",
                "Eff_allele",
                "Eff_Allele_Freq",
                "P-Value",
                "LM_PValue"
            ]
            for var in ["Intercept","Geno"] + self.covar_labels:
                for t in ["mean", "mean_stder", "mean_pval", "var", "var_stder", "var_pval"]:
                    header.append("%s_%s" % (var.lower(), t))
            print("\t".join(header), file=f)
        else:
            print("\t".join([
                "Chr",
                "Pos",
                "RSID",
                "Phenotype",
                "N",
                "Ref_allele",
                "Eff_allele",
                "Eff_Allele_Freq",
                "P-Value",
                "geno_mean",
                "geno_mean_stderr",
                "geno_mean_pval",
                "geno_var",
                "geno_var_stderr",
                "geno_var_pval",
            ]), file=f)


    def stringify(self, value):
        try:
            return "%.3e" % (value)
        except:
            return value

    def print_result(self, f=sys.stdout, verbose=False):
        """Print result to f

        :param f: stream to print output
        :param verbose: print all data or only the most important parts?
        """
        var_count = len(self.betas)/2
        if verbose:
            results = [str(x) for x in [
                self.chr,
                self.pos,
                self.rsid,
                self.ph_label,
                self.non_miss,
                self.maj_allele,
                self.min_allele,
                "%.4e" % self.eff_alcount,
                "%.4e" % self.p_mvtest,
                "%.4e" % self.lmpv
            ]]

            for i in range(var_count):
                results.append("%.3e" % self.betas[i])
                results.append(self.stringify(self.beta_stderr[i]))
                results.append(self.stringify(self.beta_pvalues[i]))
                results.append(self.stringify(self.betas[i+var_count]))
                results.append(self.stringify(self.beta_stderr[i+var_count]))
                results.append(self.stringify(self.beta_pvalues[i+var_count]))
            print("\t".join([str(x) for x in results]), file=f)

        else:
            print("\t".join(str(x) for x in [
                self.chr,
                self.pos,
                self.rsid,
                self.ph_label,
                self.non_miss,
                self.maj_allele,
                self.min_allele,
                "%.3e" % self.eff_alcount,
                "%.3e" % self.p_mvtest,
                "%.3e" % self.betas[1],
                self.stringify(self.beta_stderr[1]),
                "%.3e" % self.beta_pvalues[1],
                "%.3e" % self.betas[1+var_count],
                self.stringify(self.beta_stderr[1+var_count]),
                "%.3e" % self.beta_pvalues[1+var_count],
            ]), file=f)


