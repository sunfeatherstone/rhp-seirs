

from pathlib import Path
import pandas as pd


SRC_DIR   = Path("rhythm/raw_posts_idp")
PART_FILES = [
    SRC_DIR / "posts_idp_part1.parquet",
    SRC_DIR / "posts_idp_part2.parquet",
    SRC_DIR / "posts_idp_part3.parquet",
]
DST_FILE  = "rhythm/raw_posts_idp/posts_idp.parquet"

dfs = []
for p in PART_FILES:
    print(f"Loading {p} …")
    dfs.append(pd.read_parquet(p))

print("Concatenating …")
merged = pd.concat(dfs, ignore_index=True)

merged.to_parquet(DST_FILE, index=False)
print(f"Merged file written to {DST_FILE}")

# import pandas as pd, pathlib as pl
# orig  = pd.read_parquet("/Users/sunyushi/Downloads/public/rhythm/posts_idp.parquet")
# merged = pd.read_parquet(pl.Path("rhythm/posts_idp_merged.parquet"))
# assert len(orig) == len(merged) and orig.equals(merged)