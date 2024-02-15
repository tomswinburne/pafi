import glob
import numpy as np
from pafi import ResultsProcessor

"""
    All test data here is for the same system, currently at zero temperature
    TODO: have finite temperature tests..
"""

csv_list = glob.glob("dumps/pafi_data_*.csv")
p = ResultsProcessor(data_path=csv_list)

# Returns a pandas DataFrame and a list of dictionaries
x_key = 'ReactionCoordinate'
y_key = 'FreeEnergyGradient'
_ , integrated_data = p.integrate(argument=x_key,
                target=y_key,remesh=10,
                return_remeshed_array=True)
for int_dat in integrated_data:
    ii = int_dat[y_key+"_integrated"].argmax()
    f_b = np.round(int_dat[y_key+"_integrated"][ii],4)
    f_e = np.round(int_dat[y_key+"_integrated_err"][ii],4)    
    print(f"""
        Temperature: {int_dat['Temperature']}
        Free Energy Barrier: {f_b} +- {f_e}eV
        """)
