from pathlib import Path
import pandas as pd
SRC_FILE   = Path("rhythm/raw_posts_idp/posts_idp.parquet")
DST_DIR    = Path("rhythm/raw_posts_idp")          
BASE_NAME  = "posts_idp_part"       
N_PARTS    = 3                      

print(f"Loading {SRC_FILE} …")
df = pd.read_parquet(SRC_FILE)      

n_total = len(df)
chunk_sz = (n_total + N_PARTS - 1) // N_PARTS  
print(f"Total rows: {n_total}  ->  {N_PARTS} parts, ~{chunk_sz} rows each")

for i in range(N_PARTS):
    start = i * chunk_sz
    stop  = min((i + 1) * chunk_sz, n_total)
    if start >= stop:
        break                              
    part_df = df.iloc[start:stop].reset_index(drop=True)
    
    dst_path = DST_DIR / f"{BASE_NAME}{i+1}.parquet"
    part_df.to_parquet(dst_path, index=False)
    print(f"  Wrote rows {start}:{stop}  ->  {dst_path}")

print("✅ Split done.")