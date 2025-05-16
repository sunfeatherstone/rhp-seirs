import pandas as pd

input_parquet = "rhythm/preprocessed_data/posts_idp_interpolation.parquet"
df_in = pd.read_parquet(input_parquet)
print(df_in)
df_in['ms_time'] = pd.to_datetime(df_in['ms_time'], format="%Y-%m-%d %H:%M:%S.%f", errors='coerce')
df_in['weekday'] = df_in['ms_time'].dt.day_name()
df_in['date'] = df_in['ms_time'].dt.date
valid_dates = []
unique_dates_by_weekday = (
    df_in.groupby('weekday')['date']
      .nunique()
      .reindex(['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'])
)
print(unique_dates_by_weekday)
for i in [90,95,98,99,100]:
  df=df_in
  for day, group in df.groupby('weekday'):
      date_counts = group.groupby('date').size()
      lower_bound = date_counts.quantile(1-i/100)
      upper_bound = date_counts.quantile(i/100)
      valid = date_counts[(date_counts >= lower_bound) & (date_counts <= upper_bound)].index
      valid_dates.extend([(day, d) for d in valid])
  valid_keys = set(valid_dates)
  df = df[df.apply(lambda row: (row['weekday'], row['date']) in valid_keys, axis=1)]
  unique_dates_by_weekday = (
      df.groupby('weekday')['date']
        .nunique()
        .reindex(['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'])
  )
  output_parquet = f"rhythm/preprocessed_data/posts_idp_interpolation_truncated_{i}.parquet"
  df.to_parquet(output_parquet, index=False)
  in_path  = f"rhythm/preprocessed_data/posts_idp_interpolation_truncated_{i}.parquet"
  out_path = f"rhythm/preprocessed_data/posts_idp_interpolation_truncated_{i}_aligned.parquet"
  df = pd.read_parquet(in_path)
  print(df)
  days_cnt = df.groupby('weekday')['date'].nunique()
  min_days = days_cnt.min()
  keep_p = (min_days / days_cnt).to_dict()       # {weekday: keep_fraction}
  def downsample(g):
      p = keep_p[g.name]        # g.name == 当前 weekday
      return g if p >= 1 else g.sample(frac=p, random_state=0)
  aligned_df = (
      df.groupby('weekday', group_keys=False)
        .apply(downsample)
  )
  aligned_df.to_parquet(out_path, index=False)
  # ── freq2prob.py ─────────────────────────────────────────────────────────
  import pandas as pd, numpy as np
  from scipy.interpolate import CubicSpline
  from pathlib import Path
  from tqdm import tqdm
  import matplotlib.pyplot as plt

  plt.rcParams['font.sans-serif'] = ["PingFang HK"]
  plt.rcParams['axes.unicode_minus'] = False
  PARQUET = Path(f"rhythm/preprocessed_data/posts_idp_interpolation_truncated_{i}_aligned.parquet")
  print("Loading ms_time …")
  df = pd.read_parquet(PARQUET, columns=["ms_time"])
  df["ms_time"] = pd.to_datetime(df["ms_time"])
  df["bucket_of_week"] = df["ms_time"].dt.dayofweek * 12 + (df["ms_time"].dt.hour // 2)
  counts = (df.groupby("bucket_of_week")
              .size()
              .reindex(range(84), fill_value=0)
              .astype("int64"))
  prob = counts / counts.sum()          # ρ₂h(b)  –  shape (84,)
  hours = np.arange(85)
  prob_arr = np.append(prob.values, prob.values[0])
  cs = CubicSpline(hours, prob_arr, bc_type="periodic")
  import scipy.io as sio
  breaks = cs.x                    # shape: (85,)   0,1,…,84
  coefs  = cs.c.T                  # SciPy (4,n) → (n,4)
  pp_mat = {
      "breaks": breaks,            # 1×(n+1)
      "coefs" : coefs,             # n×4  (cubic)
      "order" : np.array([4])     
  }
  sio.savemat(f"rhythm/preprocessed_data/weibo_spline_pp_{i}.mat", pp_mat)
  counts = (
      df.groupby("bucket_of_week")
        .size()
        .reindex(range(84), fill_value=0)
        .astype("int64")
  )
  counts.to_csv(f"rhythm/preprocessed_data/bucket_counts_in_{i}.csv",
                index_label="bucket_of_week",
                header=["count"])