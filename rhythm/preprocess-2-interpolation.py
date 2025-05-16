import numpy as np
import pandas as pd
from pathlib import Path

INPUT_PATH  = Path('rhythm/raw_posts_idp/posts_idp.parquet')
OUTPUT_PATH = Path('rhythm/preprocessed_data/posts_idp_interpolation.parquet')
df = pd.read_parquet(INPUT_PATH)
df['ms_time'] = pd.to_datetime(df['ms_time'])
df['hour_start'] = pd.to_datetime(df['time_point'], format='%Y-%m-%d-%H')
np.random.seed(0)                
rows_to_add = []
grp_cols = ['keyword', 'time_point']
for (kw, tp), g in df.groupby(grp_cols, sort=False):
    if g['page'].max() != 45:
        continue                 
    hour_start  = g['hour_start'].iloc[0]  
    first_time  = g['ms_time'].min()       
    if first_time <= hour_start:
        continue           
    later       = g[g['ms_time'] > first_time]
    if len(later) >= 2:
        duration = (later['ms_time'].max() - first_time).total_seconds()
        lam      = max((len(later) - 1) / duration, 1/3600) 
    else:
        lam      = 1/60.0  
    t = first_time
    while True:
        delta  = np.random.exponential(scale=1/lam)
        t     -= pd.Timedelta(seconds=delta)
        if t <= hour_start:
            break
        rows_to_add.append({
            'keyword'      : kw,
            'time_point'   : tp,
            'page'         : 46,
            'ms_time'      : t,
            'weibo_id_rand': np.nan,
            'hour_start'   : hour_start
        })
if rows_to_add:
    df_extra = pd.DataFrame(rows_to_add)
    df       = pd.concat([df, df_extra], ignore_index=True)

df = df.sort_values(['keyword', 'time_point', 'ms_time']).reset_index(drop=True)
df = df.drop(columns=['hour_start'])
df.to_parquet(OUTPUT_PATH, index=False)
print(df)