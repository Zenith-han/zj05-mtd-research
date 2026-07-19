# 来源说明

分析使用的处理后指标文件为：

```text
data/environment_proxy_metrics.csv
```

用于构建或解释这些指标的古环境记录包括：

- 更新世相对海平面堆栈：Spratt and Lisiecki (2016)。
- ODP Site 1146 浮游有孔虫氧同位素记录：Caballero-Gill et al. (2012)。
- 中国黄土粒度堆栈：Yang and Ding (2014)。
- 中国黄土磁化率堆栈：Sun et al. (2006)。
- 中国洞穴石笋氧同位素堆栈：Cheng et al. (2016)。

指标文件包含以下核心字段：

```text
age_ka,RSL_State_Value,EASM_Intensity_Value,RSL_RateofChange_dka_h1.2,EASM_RateofChange_dka_h1.2
```

五个字段均为必需字段。两个变化率列是按正文 Equation (1) 计算的处理后指标：

```text
[X(t - 0.6 ka) - X(t + 0.6 ka)] / 1.2 ka
```

`RSL_RateofChange_dka_h1.2` 的单位为 m/ka，
`EASM_RateofChange_dka_h1.2` 的单位为 EASM 指数单位/ka。
`04_build_environment_timeseries.R` 直接读取 CSV 中的这两列；任一列缺失即中止，公开流程不重新计算 rate。

`04_build_environment_timeseries.R` 根据连续状态/强度指标，使用 `code/analysis_config.R` 定义的对称 51 点窗口计算两个变率序列。输入文件有意不保存预先计算的变率列。随后，脚本根据连续指标重建全部分类标签。
