import json
import pysam
import operator


def parse_condition(cond_str):
    for op_str in sorted(ops.keys(), key=len, reverse=True):
        if cond_str.startswith(op_str):
            value = cond_str[len(op_str):].strip()
            if not value:
                return None, None
            return ops[op_str], value
    return None, None


ops = {
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne
}

# 篩選判斷    
def variant_passes(record, criteria_list):
    for criteria in criteria_list:
        field = criteria[1]
        op_func = criteria[2]
        value = criteria[3]
        
        if criteria[0]=="INFO":            
            vcf_val = record.info.get(field)
            if vcf_val is None:
                return False
            if isinstance(vcf_val, (list, tuple)):
                vcf_val = vcf_val[0]
            try:
                if not op_func(float(vcf_val), float(value)):
                    return False
            except Exception:
                return False
                
        if criteria[0]=="FILTER":
            
            filt_val = next(iter(record.filter.keys()), "PASS")
            if not op_func(str(filt_val), str(value)):
                return False
    
    return True

# 驗證 JSON 條件結構
def validate_criteria_json(criteria):
    for section, fields in criteria.items():
        if not isinstance(fields, dict):
            raise ValueError(f"Section '{section}' should be a dictionary.")
        for key, cond_str in fields.items():
            if not isinstance(cond_str, str):
                raise ValueError(f"Condition for '{section}.{key}' must be a string.")

            op_func, value = parse_condition(cond_str)
            if op_func is None or value is None:
                raise ValueError(f"Invalid condition in '{section}.{key}': '{cond_str}'")

    print("JSON format is valid.")
    



if __name__ == "__main__":

    with open("criteria.json", "r") as f:
        criteria_json = json.load(f)
        validate_criteria_json(criteria_json)
       
    
    criteria_list=[]
    for section in criteria_json:
        for category, values in criteria_json[section].items():
            op_func, raw_val = parse_condition(values)
            try:
                val = float(raw_val)
            except ValueError:
                val = raw_val.strip('"')
            criteria_list.append([section,category,op_func,val])
            #print(criteria_list)
            
    vcf = pysam.VariantFile("SAMPLE_mutect2_raw.vcf", "r")
    vcf_out = pysam.VariantFile("SAMPLE_mutect2.vcf", "w", header=vcf.header)
    
    for record in vcf.fetch():
        if variant_passes(record, criteria_list):
            record.filter.clear()
            record.filter.add("PASS")
        vcf_out.write(record)
    
