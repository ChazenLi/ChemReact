# schemas/ - JSON Schema 定义

输入输出数据格式验证。示例数据在 `samples/`。

## 文件

| 文件 | 说明 |
|------|------|
| `host_request.schema.json` | 主机请求 schema (target_smiles 必需, 规划器/审计选项) |
| `routes.schema.json` | 路线数据 schema (route_id, steps, audit_verdict, score) |
| `strategy.schema.json` | 策略分析 schema (core_skeleton, key_disconnections, risk_budget) |
| `vis_plan.schema.json` | 可视化计划 schema (overview_grid, tree_view, step_views) |
