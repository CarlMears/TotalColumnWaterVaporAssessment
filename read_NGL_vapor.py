import numpy as np
import pandas as pd
from pathlib import Path

data_path = Path('M:/GPS_PW_Compare/python/TotalColumnWaterVaporAssessment/gnss_data/NGL_repro3_v1.0-20231013T183913Z-001/NGL_repro3_v1.0/daily/cumul_1994_2022/GPS_NGL_repro3_v1.0_cumul_1994_2022_daily/GPS')
def read_NGL_vapor(station_id: str,
                   year_to_keep: int = None):
        
    station_file = data_path / f'{station_id}.iwv'
    df = pd.read_csv(station_file,delim_whitespace=True)

    #remove the commas from the column names
    for key in df.keys():
        if ',' in key:
            new_key = key.translate({ord(i): None for i in ','})
            df.rename(columns={key:new_key},inplace=True)

    
    if year_to_keep is not None:
        df = df.loc[(df['yyyy'] == year_to_keep) & (np.isfinite(df['iwv']))]
    else:
        df = df.loc[np.isfinite(df['iwv'])]

    df.reset_index(inplace=True,drop=True)
    return df


if __name__ == '__main__':
    df = read_NGL_vapor('j521',year_to_keep=2019)