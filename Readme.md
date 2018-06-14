#### Usage:

```bash
bash clear_result.sh
bash create_r_script.sh
bash runningSessions.sh
# check session1
cat session1/test.Rout
# check all sessions
bash collect_result.sh && cat clean_data.r.Rout
```

#### Notes:

When changing parameters, you need to update the following:  

- `par, ub, lb` in `test.R`
- `ub, lb` in `conv_plot.R`
- `par` in `clean_data.R`

