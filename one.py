import pandas as pd
import math
import numpy as np
import re

def split_item(input_string):
    match = re.match(r"([A-Z][a-z]*)(\d+)", input_string)
    if match:
        element = match.group(1)
        number = int(match.group(2))
        return element, number
    else:
        return None, None


def Mixentropy(Cp):
    Sum = 0
    for i in Cp:
        Sum = Sum + float(Cp[i]) * np.log(float(Cp[i]))
    DeltaS = -8.3144621 * Sum
    return DeltaS

def ElecDiff(Cp , element_data):
    Sum = 0
    xbar = 0
    for i in Cp:
        xbar = xbar + float(Cp[i]) * float(element_data.loc[element_data['Symbol'] == i, 'Electronegativity'].item())
    for i in Cp:
        Sum = Sum + float(Cp[i]) * np.power((float(element_data.loc[element_data['Symbol'] == i, 'Electronegativity'].item()) - xbar), 2)
    ED = np.sqrt(Sum)
    return ED
def FVEC(Cp , element_data):
    VEC = 0
    for i in Cp:
        VEC = VEC + float(Cp[i]) * float(element_data.loc[element_data['Symbol'] == i, 'VEC'].item())
    return VEC
def normalize_composition(composition):
    total = sum(composition.values())
    return {element: amount / total for element, amount in composition.items()}

def get_enthalpy(E1, E2):
    try:
        Em = enthalpy_data.loc[enthalpy_data['Symbol'] == E1, E2].values[0]
        if np.isnan(Em):  # If NaN, check the reverse pair
            Em = enthalpy_data.loc[enthalpy_data['Symbol'] == E2, E1].values[0]
        return Em
    except (IndexError, KeyError):
        raise ValueError(f"No enthalpy data found for pair ({E1}, {E2})")

def calculate_mixing_enthalpy(normalized_composition):
    elements = list(normalized_composition.keys())
    Emix = 0
    for i in range(len(elements) - 1):
        for j in range(i + 1, len(elements)):
            E1, E2 = elements[i], elements[j]
            Em = get_enthalpy(E1, E2)
            Emix += 4 * Em * normalized_composition[E1] * normalized_composition[E2]
    return Emix
def AtmSizeDiff(Cp , element_data):
    Sum = 0
    rbar = 0
    for i in Cp:
        rbar = rbar + float(Cp[i]) * float(element_data.loc[element_data['Symbol'] == i, 'AtomicRadius'].item())
    for i in Cp:
        Sum = Sum + float(Cp[i]) * np.power(
            (1 - float(element_data.loc[element_data['Symbol'] == i, 'AtomicRadius'].item()) / rbar), 2)
    ASD = 100 * np.sqrt(Sum)
    return ASD
def calculate_melting_points(data, elements_dict):

    melting_point_vals = []

    for element_name, composition in elements_dict.items():
        # composition /= 100  
        
        filtered_data = data[data["Symbol"].str.lower() == element_name.lower()]
        
        if not filtered_data.empty:
            d = filtered_data.iloc[0].to_dict()
            melting_point_key = [key for key in d.keys() if "MeltingPoint" in key]

            if melting_point_key:
                m_point = d[melting_point_key[0]]
                melting_point_vals.append(m_point)
            else:
                print(f"Melting point not found for {element_name}")
        else:
            print(f"Element '{element_name}' not found in the dataset.")

    if melting_point_vals:
        final_meltingpoint = sum(melting_point_vals) / len(melting_point_vals)
        lower, upper = 0.6 * final_meltingpoint, 0.8 * final_meltingpoint

        return final_meltingpoint, lower, upper
    else:
        print("No melting points were recorded.")
        return 0


def calculate_md_values(elements_dict, data):

    md_values_df = data.copy()
    total_md_value = 0

    for element_name, composition in elements_dict.items():
        # composition /= 100

        filtered_data = md_values_df[md_values_df["Element"].str.lower() == element_name.lower()]

        if not filtered_data.empty:
            md_val = filtered_data.iloc[0]["MD_value"]
            total_md_value += composition * md_val
        else:
            print(f"MD value not found for {element_name}")

    crystal = None
    if total_md_value < 0.891:
        crystal = "FCC"
    elif total_md_value >= 0.891:
        crystal = "Intermetallic + FCC"
        
    return crystal, total_md_value

def parse_elements(s):
    elements = s.split(',')
    ele_dict = {}
    
    for element in elements:
        symbol = ''.join([ch for ch in element if ch.isalpha()])  # Extract the symbol
        value = int(''.join([ch for ch in element if ch.isdigit()]))  # Extract the numeric value
        ele_dict[symbol] = value / 100  
    
    return ele_dict


print("Enter elements separated by commas ")
print("For example:  -->  Fe40,Ni50,Cr10  ")

s=input()

ele = parse_elements(s)

print(ele)

elementdata    = pd.read_csv("elementdata.csv")
enthalpy_data  = pd.read_csv("enthalpydata.csv")  
melting_points = pd.read_csv("meltingpoint.csv")
md_values      = pd.read_csv("md_vals.csv")


res  = AtmSizeDiff(ele, elementdata) 
res1 = FVEC(ele, elementdata) 
res2 = ElecDiff(ele, elementdata) 
res3 = Mixentropy(ele)

final_meltingpoint, lower, upper = calculate_melting_points( melting_points , ele )
crystal , mdval  = calculate_md_values( ele , md_values)



normalized = normalize_composition(ele)
mixing_enthalpy = calculate_mixing_enthalpy(normalized)


print(f"Atomic size difference: {res:.4f}")
print(f"Electronegativity difference: {res2:.4f}")
print(f"(VEC) Valence electron concentration : {res1:.4f}")

print()

print(f"(ΔSmix) Mixing Entropy : {res3:.4f}")
print(f'(ΔHmix) Mixing Enthalpy : {mixing_enthalpy:.3f}')

print(f"Crystal : {crystal}")
print(f"MD-Value : {mdval} ")

print()

print(f"Average melting point: {final_meltingpoint}")
print(f"Lower limit: {lower}")
print(f"Upper limit: {upper}")