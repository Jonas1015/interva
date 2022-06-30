# -*- coding: utf-8 -*-

"""
interva.interva5
-------------------

This module contains the class for the InterVA5 algorithm.
"""

from pandas import (DataFrame, Series, read_csv, read_excel, to_numeric, isna)
from numpy import (full, ndarray, nan, nansum, array, amax, argmax, delete, 
                   hstack, savetxt, zeros)
from os import path, chdir, getcwd, mkdir
from logging import basicConfig, info, warning
from csv import writer
from pkgutil import get_data
from io import BytesIO

from interva.data.causetext import CAUSETEXTV5
from vacheck.datacheck5 import datacheck5


class InterVA5:
    """InterVA5 algorithm for assigning cause of death.

    :param va_input: Verbal Autopsy data
    :type va_input: pandas data.frame or path to CSV file
    :param hiv: likelihood of HIV as a cause of death.  Possible values are
     "H" for high (~ 1:100 deaths), "L" for low (~ 1:1000), or "V" for very
     low (~ 1:10000)
    :type hiv: string
    :param malaria: likelihood of malaria as a cause of death.  Possible values are
     "H" for high (~ 1:100 deaths), "L" for low (~ 1:1000), or "V" for very
     low (~ 1:10000)
    :type malaria: string
    :param write: logical value
    :type write: boolean
    :param directory: The directory to store the output from InterVA5.
     It should either be an existing valid directory, or a new folder to be created.
     If no path is given and the parameter for "write" is True, then
     the function stops and an error message is produced.
    :type directory: directory or string
    :param filename: the filename the user wishes to save the output.
     No extension needed. The output is in .csv format by default.
    :type filename: string
    :param output: the format of the output. Possible Values are 
     "classic": the same deliminated output format as InterVA5, or
     "extended": delimited output followed by full distribution of cause of death
     probability
    :type output: string
    :param append: a logical value indicating whether or not the new output should
     be appended to the existing file.
    :type append: boolean
    :param groupcode: a logical value indicating whether or not the group code will
     be included in the output causes.
    :type groupcode: boolean
    :param sci: an array containing the symptom-cause-information (aka Probbase)
     that InterVA uses to assign a cause of death.
    :type sci: pandas data.frame or numpy ndarray
    :param return_checked_data: a logical value indicating if the checked data
     (i.e., the data that have been modified by the consistency checks) should
     be returned.
    :type return_checked_data: boolean
    :param others: not used
    """

    def __init__(self, va_input, hiv: str, malaria: str, write: bool = True, 
                 directory = None, filename: str = "VA5_result", 
                 output: str = "classic", append: bool = False, 
                 groupcode: bool = False, sci = None, 
                 return_checked_data: bool = False, *others) -> dict:

        self.va_input = va_input
        self.hiv = hiv
        self.malaria = malaria
        self.write = write
        self.directory = directory
        self.filename = filename
        self.output = output
        self.append = append
        self.groupcode = groupcode
        self.sci = sci
        self.return_checked_data = return_checked_data
        basicConfig(filename="errorlogV5.txt", filemode='a')

    def _check_data(self, va_input: Series, va_id: str, 
                    insilico_check: bool = False):
        """Run data check."""
        
        return datacheck5(va_input, va_id, insilico_check)
        

    def run(self):
        """Assign causes of death.
        
        :return: ids of VA input, VA results with cause assignments and
         likelihoods, likelihood of malaria and HIV as causes of death,
         and cleaned data from data consistency checks.
        :rtype: dictionary with keys ID (pandas.series), 
         VA_result (pandas data.frame), Malaria (str), HIV (str), and
         checked_data (pandas data.frame).
        """
        
        def va5(id, malprev, hivprev, pregstat, preglik, cause1, lik1, cause2, 
                lik2, cause3, lik3, indet, comcat, comnum, wholeprob, *others):
            return [id, str(malprev), str(hivprev), pregstat, preglik, 
                    cause1, lik1, cause2, lik2, cause3, lik3, indet, 
                    str(comcat), comnum, wholeprob]
                    
        def save_va5(x: ndarray, filename: str, write: bool):
            if not write:
                return()
            x = delete(x, 14)
            filename = filename + '.csv'
            with open(filename, 'a') as csvfile:
                savetxt(csvfile, x, delimiter=",", fmt="%f", header="")
        
        def save_va5_prob(x: ndarray, filename: str, write: bool):
            if not write:
                return()
            prob = x[14]
            x = delete(x, 14)
            filename = filename + ".csv"
            x = hstack(x, prob).squeeze()
            with open(filename, 'a') as csvfile:
                savetxt(csvfile, x, delimiter=",", fmt="%f", header="")
        
        if self.directory is None and self.write:
            raise IOError \
                ("error: please provide a directory (required when write = True)")
        if self.directory is None:
            self.directory = getcwd()
        if not path.isdir(self.directory):
            mkdir(self.directory)
        globle_dir = getcwd()
        chdir(self.directory)
        
        probbaseV5 = None
        if self.sci is None:
            probbase_xls = get_data("interva", "data/probbase.xls")
            probbase = read_excel(probbase_xls)
            probbaseV5 = probbase.to_numpy()
            probbaseV5 = delete(probbaseV5, 0, axis=0)
            self.probbaseV5Version = probbaseV5[0, 2]
        if self.sci is not None:
            validSCI = True
            if not isinstance(self.sci, DataFrame) and \
                not isinstance(self.sci, ndarray):
                validSCI = False
            if self.sci.shape[0] != 354 or self.sci.shape[1] != 87:
                validSCI = False
            if not validSCI:
                raise IOError \
                    ("Error: Invalid SCI (must be Pandas DataFrame or \
                     Numpy ndarray with 354 rows and 87 columns).")
            if isinstance(self.sci, DataFrame):
                self.sci = self.sci.to_numpy()
            probbaseV5 = self.sci
            self.probbaseV5Version = probbaseV5[0, 2]
        print("Using Probbase version:", self.probbaseV5Version, sep=" ")
        causetextV5_horizontal = DataFrame(CAUSETEXTV5)
        self.causetextV5 = causetextV5_horizontal.transpose()
        if self.groupcode:
            self.causetextV5.drop(self.causetextV5.columns[0], axis=1)
        else:
            self.causetextV5.drop(self.causetextV5.columns[1], axis=1)
        if self.write:
            info("Error & warning log built for InterVA5")
        if isinstance(self.va_input, str) and self.va_input[-4:] == ".csv":
            self.va_input = read_csv(self.va_input)
        if "i183o" in self.va_input.columns:
            self.va_input.rename(columns={'i183o': 'i183a'}, axis='columns', 
                                 inplace=True)
            print("Due to the inconsistent names in the early version of " + 
                  "InterVA5, the indicator 'i183o' has been renamed as 'i183a'.")
        
        va_data = self.va_input
        va_input_names = va_data.columns
        id_inputs = va_data.iloc[:, 0]
        va_data = va_data.to_numpy()
        if va_data.shape[0] < 1:
            raise IOError("error: no data input")
        N = va_data.shape[0]
        S = va_data.shape[1]
        if S != probbaseV5.shape[0]:
            raise IOError \
                ("error: invalid data input format. Number of values incorrect")
        if va_input_names[S-1].lower() != "i459o":
            raise IOError("error: the last variable should be 'i459o'")
        va_data_csv = get_data("interva", "data/randomva5.csv")
        randomVA5 = read_csv(BytesIO(va_data_csv))
        valabels = randomVA5.columns
        count_changelabel = 0
        for i in range(S):
            input_col = va_input_names[i]
            std_col = valabels[i]
            if input_col.lower() != std_col.lower():
                warning("Input column '{input_col}' does not match \
                        InterVA5 standard: '{std_col}'")
                count_changelabel = count_changelabel + 1
        if count_changelabel > 0:
            warning("{count_changelabel} column names changed in input.\n" + 
                    "If the change is undesirable, please change in the input " +
                    "to match standard InterVA5 input format.")
            va_input_names = valabels
        prob_ncols = probbaseV5.shape[1]
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "I"] = 1
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "A+"] = 0.8
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "A"] = 0.5
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "A-"] = 0.2
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "B+"] = 0.1
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "B"] = 0.05
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "B-"] = 0.02
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "C+"] = 0.01
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "C"] = 0.005
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "C-"] = 0.002
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "D+"] = 0.001
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "D"] = 5e-04
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "D-"] = 1e-04
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "E"] = 1e-05
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == "N"] = 0
        probbaseV5[:,17:prob_ncols][probbaseV5[:,17:prob_ncols] == ""] = 0
        probbaseV5[0, 0:17] = 0
        Sys_Prior = to_numeric(probbaseV5[0, :])
        D = len(Sys_Prior)
        self.hiv = self.hiv.lower()
        self.malaria = self.malaria.lower()
        if self.hiv not in ['h', 'l', 'v'] or self.malaria not in ['h', 'l', 'v']:
            raise IOError("error: the HIV and Malaria indicator " +
                          "should be one of the three: 'h', 'l', 'v'")
        if self.hiv == "h":
            Sys_Prior[22] = 0.05
        if self.hiv == "l":
            Sys_Prior[22] = 0.005
        if self.hiv == "v":
            Sys_Prior[22] = 1e-05
        if self.malaria == "h":
            Sys_Prior[24] = 0.05
            Sys_Prior[44] = 0.05
        if self.malaria == "l":
            Sys_Prior[24] = 0.005
            Sys_Prior[44] = 1e-05
        if self.malaria == "v":
            Sys_Prior[24] = 1e-05
            Sys_Prior[44] = 1e-05
            
        ID_list = []
        VA_result = [[] for _ in range(N)]
        if self.write and not self.append:
            header = ["ID", "MALPREV", "HIVPREV", "PREGSTAT", "PREGLIK", 
                      "CAUSE1", "LIK1", "CAUSE2", "LIK2", "CAUSE3", "LIK3", 
                      "INDET", "COMCAT", "COMNUM"]
            if self.output == "extended":
                header = header.append(self.causetextV5.iloc[:, 1])
                with open(self.filename + ".csv", 'a+', newline='') as write_obj:
                    csv_writer = writer(write_obj)
                    csv_writer.writerow(header)
        nd = max(1, round(N/100))
        np = max(1, round(N/10))
        
        if self.write:  
            info("\n\nthe following records are incomplete and " +
                "excluded from further processing:\n\n")
            
        first_pass = []
        second_pass = []
        errors = None
        if self.return_checked_data:
            self.checked_data = [[] for _ in range(N)]
        for i in range(N):
            k = i + 1
            if k % nd == 0:
                print(".", end="")
            if k % np == 0:
                percentage = round(k/N * 100)
                print(percentage, "% completed", sep="")
            if k == N:
                print("100% completed")
            index_current = str(va_data[i, 0])
            index_current = id_inputs[i]
            va_data[i, :][va_data[i, :] == "n"] = "0"
            va_data[i, :][va_data[i, :] == "N"] = "0"
            va_data[i, :][va_data[i, :] == "y"] = "1"
            va_data[i, :][va_data[i, :] == "Y"] = "1"
            for j in range(va_data.shape[1]):
                if va_data[i, j] != "0" and va_data[i, j] != "1":
                    va_data[i, j] = nan
            input_current = va_data[i, :]
            input_current[:][input_current[:] == "0"] = 0
            input_current[:][input_current[:] == "1"] = 1
            to_numeric(input_current[:])
            
            input_current[0] = 0
            if nansum(input_current[5:12]) < 1:
                if self.write:
                    errors = (errors + index_current + 
                              " Error in age indicator: Not Specified")
                continue
            if nansum(input_current[3:5]) < 1:
                if self.write:
                    errors = (errors + index_current + 
                              " Error in sex indicator: Not Specified")
                continue
            if nansum(input_current[20:328]) < 1:
                if self.write:
                    errors = (errors + index_current + 
                              " Error in indicators: No symptoms specified")
                continue
            
            input_current = Series(input_current, index=va_input_names)
            tmp = datacheck5(va_input=input_current, va_id=index_current)
            if self.return_checked_data:
                current_checked = [[id_inputs[i]]]
                current_checked.append(list(tmp["output"][1:S]))
                current_checked = \
                    [item for sublist in current_checked for item in sublist]
                self.checked_data[i] = current_checked
            input_current = tmp["output"]
            first_pass.append(tmp["first_pass"])
            second_pass.append(tmp["second_pass"])
            
            subst_vector = full(S, nan)
            subst_vector[probbaseV5[:, 5] == "N"] = 0
            subst_vector[probbaseV5[:, 5] == "Y"] = 1
            
            new_input = zeros(S, dtype=int)
            for y in range(1,S):
                if not isna(input_current[y]):
                    if input_current[y] == subst_vector[y]:
                        new_input[y] = 1
            
            input_current[input_current == 0] = 1
            input_current[0] = 0
            input_current[isna(input_current)] = 0
            reproductiveAge = 0
            preg_state = " "
            lik_preg = " "
            if input_current[4] == 1 and \
                (input_current[16] == 1 or input_current[17:19].any() == 1):
                reproductiveAge = 1
            prob = Sys_Prior[17:D]
            input_len = len(input_current)
            temp = new_input[1:input_len][new_input[1:input_len] == 1]
            for temp_sub in temp:
                for j in range(17, D):
                    prob[j-17] = prob[j-17] * \
                        to_numeric(probbaseV5[temp_sub + 1, j])
                temp_sum = sum(prob[0:3])
                temp_sum_ndarray = array(temp_sum)
                if temp_sum > 0:
                    prob[0:3] = prob[0:3] / temp_sum_ndarray
                temp_sum = sum(prob[3:64])
                temp_sum_ndarray = array(temp_sum)
                if temp_sum > 0:
                    prob[3:64] = prob[3:64] / temp_sum_ndarray
                temp_sum = sum(prob[64:70])
                temp_sum_ndarray = array(temp_sum)   
                if temp_sum > 0:
                    prob[64:70] = prob[64:70] / temp_sum_ndarray
            prob_names = self.causetextV5.iloc[:, 1]
            prob_A = prob[0:3]
            prob_B = prob[3:64]
            prob_C = prob[64:70]
            
            # Determine Preg_State and Likelihood
            prob_A_sum = sum(prob_A)
            if prob_A_sum == 0 or reproductiveAge == 0:
                preg_state = "n/a"
                lik_preg = " "
            if max(prob_A) < 0.1 and reproductiveAge == 1:
                preg_state = "indeterminate"
                lik_preg = " "
            prob_A_max_loc = argmax(prob_A)
            if prob_A_max_loc == 0 and prob_A[0] >= 0.1 and reproductiveAge == 1:
                preg_state = "Not pregnant or recently delivered"
                lik_preg = int(round(prob_A[1] / prob_A_sum * 100))
            if prob_A_max_loc == 1 and prob_A[1] >= 0.1 and reproductiveAge == 1:
                preg_state = "Pregnancy ended within 6 weeks of death"
                lik_preg = int(round(prob_A[1] / prob_A_sum * 100))
            if prob_A_max_loc == 2 and prob_A[2] >= 0.1 and reproductiveAge == 1:
                preg_state = "Pregnant at death"
                lik_preg = int(round(prob_A[2] / prob_A_sum * 100))
            
            # Determine the output of InterVA
            prob_temp = prob_B
            prob_temp_names = prob_names[3:64]
            prob_C_names = prob_names[64:70]
            top3 = [[] for _ in range(3)]
            cause1 = lik1 = cause2 = lik2 = cause3 = lik3 = None
            indet = 0
            if amax(prob_B) < 0.4:
                cause1 = lik1 = cause2 = lik2 = cause3 = lik3 = " "
                indet = 100
            if amax(prob_temp) >= 0.4:
                max1_loc = argmax(prob_temp)
                lik1 = round(amax(prob_temp) * 100)
                cause1 = prob_temp_names[max1_loc]
                prob_temp = delete(prob_temp, max1_loc)
                prob_temp_names.drop(prob_temp_names.index[max1_loc])
                top3.append(lik1)
                
                max2_loc = argmax(prob_temp)
                lik2 = round(amax(prob_temp) * 100)
                cause2 = prob_temp_names[max2_loc]
                if amax(prob_temp) < 0.5 * amax(prob_B):
                    lik2 = cause2 = " "
                prob_temp = delete(prob_temp, max2_loc)
                prob_temp_names.drop(prob_temp_names.index[max2_loc])
                top3.append(lik2)
                
                max3_loc = argmax(prob_temp)
                lik3 = round(amax(prob_temp) * 100)
                cause3 = prob_temp_names[max3_loc]
                if amax(prob_temp) < 0.5 * amax(prob_B):
                    lik3 = cause3 = " "
                top3.append(lik3)
                top3 = [float(x) if x != " " else 0 for x in top3]
                top3 = array(top3)
                indet = round(100 - nansum(top3))
            
            # Determine the Circumstances of Mortality CATegory (COMCAT) 
            # and probability
            comcat = ""
            comnum = None
            prob_c_max = amax(prob_C)
            if sum(prob_C) > 0:
                prob_C = prob_C / sum(prob_C)
            if prob_c_max < 0.5:
                comcat = "Multiple"
                comnum = " "
            else:
                comcat = prob_C_names[argmax(prob_C)]
                comnum = round(prob_c_max * 100)
            
            ID_list.append(index_current)
            combined_prob = Series(prob, index=prob_names)
            VA_result[i] = va5(index_current, self.malaria, self.hiv, preg_state, 
                               lik_preg, cause1, lik1, cause2, lik2, cause3, lik3, 
                               indet, comcat, comnum, wholeprob=combined_prob)
            if self.output == "classic":
                save_va5(VA_result[i], filename=self.filename, write=self.write)
            if self.output == "extended":
                save_va5_prob(VA_result[i], filename=self.filename, 
                              write=self.write)
        if self.write:
            info("\n the following data discrepancies were identified and " +
                 "handled: \n" + first_pass + "\nSecond pass\n" + second_pass)
        
        chdir(globle_dir)
        if not self.return_checked_data:
            self.checked_data = "return_checked_data = False"
        else:
            self.checked_data = DataFrame(self.checked_data)
            self.checked_data.columns = va_input_names
        
        ID_list = Series(ID_list, name="ID")
        ID_nan_boolean = ID_list.isna()
        nan_indices = []
        index = 0
        for val in ID_nan_boolean:
            if val == True:
                nan_indices.append(index)
            index += 1
        ID_list = ID_list.drop(nan_indices)
        
        VA_result = DataFrame(VA_result)
        VA_result.columns = ["ID", "MALPREV", "HIVPREV", "PREGSTAT", "PREGLIK", 
                             "CAUSE1", "LIK1", "CAUSE2", "LIK2", "CAUSE3", "LIK3", 
                             "INDET", "COMCAT", "COMNUM", "WHOLEPROB"]
        VA_result = VA_result.drop(nan_indices, axis=0)
        
        self.out = {"ID": ID_list,
                    "VA5": VA_result,
                    "Malaria": self.malaria,
                    "HIV": self.hiv,
                    "checked_data": self.checked_data}
        return self.out
         
    def get_hiv(self):
        """Get HIV parameter."""
        
        print(f"HIV parameter is {self.hiv}")
        return self.hiv

    def get_malaria(self):
        """Get malaria parameter."""
        
        print(f"Malaria parameter is {self.malaria}")
        return self.malaria

    def set_hiv(self, hiv_level):
        """Set HIV parameter."""
        
        hiv_lvl = hiv_level.lower()
        if hiv_lvl in ["h", "l", "v"]:
            self.hiv = hiv_lvl
        else:
            print(f"The provided HIV level \"{hiv_level}\" is invalid.")
        return self.hiv
        print(f"HIV parameter is {self.hiv}")

    def set_malaria(self, malaria_level):
        """Set malaria parameter."""
        
        malaria_lvl = malaria_level.lower()
        if malaria_lvl in ["h", "l", "v"]:
            self.malaria = malaria_lvl
        else:
            print(f"The provided malaria level \"{malaria_level}\" is invalid.")
        return self.malaria
        print(f"Malaria parameter is {self.malaria}")

    def get_ids(self):
        """Return pandas series of ID column in data."""
        
        va_df = self.va_input
        if isinstance(va_df, str) and va_df[-4:] == ".csv":
            va_df = read_csv(va_df)
        return va_df.loc[:, "ID"]

    def plot_csmf(self, top=10, file=None):
        """Plot cause-specific mortality fraction (CSMF)."""
        pass

    def get_csmf(self, top=10):
        """Print top causes in cause-specific mortality fraction (CSMF)."""
        pass

    def write_csmf(self):
        """Write cause-specific mortality fraction (CSMF) to CSV file."""
        pass

    def get_top_cod(self, top=3, include_propensities=False):
        """Get top causes of death for each individual."""
        pass

    def write_top_cod(self, top=3, include_propensities=False):
        """Write top causes of death for each individual to CSV file."""

    def get_indiv_prob(self, include_propensities=False):
        """Get individual causes of death distribution."""
        pass

    def write_indiv_prob(self, include_propensities=False):
        """Write individual cause of death distribution to CSV file."""
        pass
