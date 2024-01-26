from interva.interva5 import InterVA5
from pycrossva.transform import transform
from interva.utils import csmf
import pandas as pd
import numpy as np

def runCCVA(odk_raw, id_col=None, instrument='2016WHOv151', algorithm='InterVA5', top=10, undetermined=True, malaria="h", hiv="h"):
    if id_col:
        input_data = transform((instrument, algorithm), odk_raw, raw_data_id=id_col)
    else:
        input_data = transform((instrument, algorithm), odk_raw)
    #input_data = transform(("2016WHOv151", "InterVA5"), odk_raw, raw_data_id="vaid")
    
    iv5out = InterVA5(input_data, hiv=hiv, malaria=malaria, write=True, directory=".", filename="ccva_results")

    iv5out.run()

    ## This is the direct method of getting total csmf all cause all groups
    # all = {
    #     "index":iv5out.get_csmf(top=top).index.tolist(),
    #     "values":iv5out.get_csmf(top=top).tolist()
    # }

    # for the below options, use the following combinations
    # age = 'adult, child, neonate'
    # sex = 'male,female'


    all = {
        "index":csmf(iv5out, top=top, age=None, sex=None).index.tolist(),
        "values":csmf(iv5out, top=top, age=None, sex=None).tolist()
    }
    if not undetermined:
        top +=1
        index = np.array(csmf(iv5out, top=top, age=None, sex=None).index.tolist())
        values = np.array(csmf(iv5out, top=top, age=None, sex=None).tolist())
        idx = np.argwhere(index == "Undetermined")
        
        if len(idx) != 0:
            index = np.delete(index,idx)
            values = np.delete(values,idx)
            
            all = {
                "index":index.tolist(),
                "values":values.tolist()
            }
    
    male = {
        "index":csmf(iv5out, top=top, age=None, sex='male').index.tolist(),
        "values":csmf(iv5out, top=top, age=None, sex='male').tolist()
    }
    if not undetermined:
        top +=1
        index = np.array(csmf(iv5out, top=top, age=None, sex='male').index.tolist())
        values = np.array(csmf(iv5out, top=top, age=None, sex='male').tolist())
        idx = np.argwhere(index == "Undetermined")
        
        
        if len(idx) != 0:
            index = np.delete(index,idx)
            values = np.delete(values,idx)
            
            male = {
                "index":index.tolist(),
                "values":values.tolist()
            }
    
    female = {
        "index":csmf(iv5out, top=top, age=None, sex='female').index.tolist(),
        "values":csmf(iv5out, top=top, age=None, sex='female').tolist()
    }
    if not undetermined:
        top +=1
        index = np.array(csmf(iv5out, top=top, age=None, sex='female').index.tolist())
        values = np.array(csmf(iv5out, top=top, age=None, sex='female').tolist())
        idx = np.argwhere(index == "Undetermined")
        
        if len(idx) != 0:
            index = np.delete(index,idx)
            values = np.delete(values,idx)
            
            female = {
                "index":index.tolist(),
                "values":values.tolist()
            }
    
    adult = {
        "index":csmf(iv5out, top=top, age='adult', sex=None).index.tolist(),
        "values":csmf(iv5out, top=top, age='adult', sex=None).tolist()
    }
    if not undetermined:
        top +=1
        index = np.array(csmf(iv5out, top=top, age='adult', sex=None).index.tolist())
        values = np.array(csmf(iv5out, top=top, age='adult', sex=None).tolist())
        idx = np.argwhere(index == "Undetermined")
        
        if len(idx) != 0:
            index = np.delete(index,idx)
            values = np.delete(values,idx)
            
            adult = {
                "index":index.tolist(),
                "values":values.tolist()
            }
    
    child = {
        "index":csmf(iv5out, top=top, age='child', sex=None).index.tolist(),
        "values":csmf(iv5out, top=top, age='child', sex=None).tolist()
    }

    if not undetermined:
        top +=1
        index = np.array(csmf(iv5out, top=top, age='child', sex=None).index.tolist())
        values = np.array(csmf(iv5out, top=top, age='child', sex=None).tolist())
        idx = np.argwhere(index == "Undetermined")
        
        if len(idx) != 0:
            index = np.delete(index,idx)
            values = np.delete(values,idx)
            
            child = {
                "index":index.tolist(),
                "values":values.tolist()
            }


    neonate = {
        "index":csmf(iv5out, top=top, age='neonate', sex=None).index.tolist(),
        "values":csmf(iv5out, top=top, age='neonate', sex=None).tolist()
    }
    if not undetermined:
        top +=1
        index = np.array(csmf(iv5out, top=top, age='neonate', sex=None).index.tolist())
        values = np.array(csmf(iv5out, top=top, age='neonate', sex=None).tolist())
        idx = np.argwhere(index == "Undetermined")
        
        if len(idx) != 0:
            index = np.delete(index,idx)
            values = np.delete(values,idx)
            
            neonate = {
                "index":index.tolist(),
                "values":values.tolist()
            }

    merged_df = pd.concat([csmf(iv5out, top=top, age=None, sex=None), 
                           csmf(iv5out, top=top, age=None, sex='male'),
                           csmf(iv5out, top=top, age=None, sex='female'),
                           csmf(iv5out, top=top, age='adult', sex=None),
                           csmf(iv5out, top=top, age='child', sex=None),
                           csmf(iv5out, top=top, age='neonate', sex=None)], axis=1)
    
    merged_df.columns = ['all','male','female','adult','child','neonate']
    merged_df.fillna(0, inplace=True)

    merged_arr = []
    for col in merged_df.columns:
        merged_arr.append(merged_df[col].tolist())

    merged = {
        "index": merged_df.index.tolist(),
        "values": merged_arr
    }

    ccva_results = {
        "all":all,
        "male":male,
        "female":female,
        "adult":adult,
        "child":child,
        "neonate":neonate,
        "merged":merged
    }

    return ccva_results


path = "WHOVA_V1_5_3_TZV1_2.csv"
df_odk = pd.read_csv(path, low_memory=False)
idColumn = 'KEY'
instrument = "2016WHOv151"
malaria = 'h'
hiv = 'h'
algorithm = 'InterVA5'
topN = 10
undetermined = True

if idColumn != "0":
    if df_odk[idColumn].squeeze().is_unique:
        output = runCCVA(df_odk, idColumn, instrument=instrument, top=topN, undetermined=undetermined, hiv=hiv, malaria=malaria)
    else:
        print("id column not unique")
else:
    output = runCCVA(df_odk, instrument=instrument, top=topN, undetermined=undetermined, hiv=hiv, malaria=malaria)

print(output)