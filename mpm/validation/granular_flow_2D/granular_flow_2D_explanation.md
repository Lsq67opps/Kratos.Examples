# granular_flow_2D 案例代码行解读（按主函数执行顺序）

本文按 `MainKratos.py` 的主函数执行顺序，对 `mpm/validation/granular_flow_2D` 中涉及的脚本与配置逐段说明。涉及计算的部分给出对应的公式，格式均使用 `[ ... ]` 便于直接渲染。

## 主入口：`source/MainKratos.py`
- **L1**：导入核心命名空间 `KratosMultiphysics`，提供模型、参数等基础类。
- **L2**：从 `KratosMultiphysics.MPMApplication` 导入 `MPMAnalysis`，这是材料点法 (MPM) 的高层求解驱动。
- **L9-L13**：在脚本作为主程序时，读取 `ProjectParameters.json` 并构造成 `Parameters` 对象。
- **L14-L16**：构造空的 `Model`，用参数实例化 `MPMAnalysis` 并调用 `simulation.Run()`。`Run()` 内部完成：读取网格/粒子模型、创建求解器、初始化过程列表、进入时间步循环、输出结果。时间积分采用 Newmark 隐式算法，其更新公式为：
  - `[ \mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t\,\mathbf{v}_n + \Delta t^2\!\left(\tfrac{1}{2}-\beta\right)\mathbf{a}_n + \beta\,\Delta t^2\,\mathbf{a}_{n+1} ]`
  - `[ \mathbf{v}_{n+1} = \mathbf{v}_n + \Delta t\!\left((1-\gamma)\mathbf{a}_n + \gamma\,\mathbf{a}_{n+1}\right) ]`
  其中 `beta, gamma` 由 Kratos 的 `newmark` 配置给出，`Δt` 为时间步长。

## 参数：`source/ProjectParameters.json`
- **L2-L8（problem_data）**：定义问题名、并行类型、时间范围 `t∈[0,2]`。给定时间步长 `Δt = 5×10^{-5}` 秒，对应步数 `[ N_{\text{step}} = \frac{2.0-0.0}{5\times 10^{-5}} ]`。
- **L10-L44（solver_settings）**：指定动态隐式求解器、2D 域、非线性迭代与收敛阈值。`model_import_settings` / `grid_model_import_settings` 分别加载材料点云 `granular_flow_2D_Body.mdpa` 与背景网格 `granular_flow_2D_Grid.mdpa`。线性求解器选择稀疏 LU。
- **L24-L26（time_stepping）**：时间步配置，对应上式中的 `Δt`。
- **L33-L35（子模型部件列表）**：`Parts_Parts_Auto1/2` 为材料点与网格的子模型部件，`DISPLACEMENT_Displacement_Auto1` 用于施加边界约束。
- **L47-L70（processes：约束与重力）**：
  - 约束：`DISPLACEMENT` 三个方向全约束 => 位移边界条件 `[ \mathbf{u} = \mathbf{0} ]`。
  - 重力：对材料点施加体加速度 `g=9.81`，方向 `(0,-1,0)`，公式 `[ \mathbf{a}_g = 9.81\,(0,-1,0) ]`，进入动量平衡 `[ \rho\,\mathbf{a} = \nabla\!\cdot\!\boldsymbol{\sigma} + \rho\,\mathbf{a}_g ]`。
- **L73-L157（output_processes）**：GiD 与 VTK 输出分别针对材料点与背景网格，按照时间步 `0.01 s` 写出。无额外计算，仅配置结果变量。

## 材料：`source/ParticleMaterials.json`
- **L1-L22**：唯一的属性集绑定到 `Initial_MPM_Material.Parts_Parts_Auto1`。选用 Hencky 模型的摩尔-库仑平面应变塑性本构。
  - 密度 `[ \rho = 2650\,\text{kg/m}^3 ]`，厚度 `[ t = 1.0 ]`，每单元材料点数 `3`。
  - 本构屈服函数（摩尔-库仑）示意：`F = 0` 满足 `[ F = \sqrt{3J_2} + p\,\tan\varphi - c ]`，其中 `[ \varphi = 0.3456\,\text{rad} ]`、`[ c = 0 ]`。
  - 弹性部分：杨氏模量与泊松比给出应力-应变关系 `[ \boldsymbol{\sigma} = \mathbf{C}:\boldsymbol{\varepsilon} ]`，`C` 为各向同性 2D 平面应变弹性矩阵。

## 材料点云：`source/granular_flow_2D_Body.mdpa`
- **L1-L6**：`ModelPartData` 与空 `Properties`，实际属性由材料文件绑定。
- **L7-L5161（Nodes）**：列举所有材料点坐标 `(id, x, y, z)`，构成初始颗粒柱。每个点的质量由 `[ m_p = \rho\,V_p ]`（体积来自生成工具）确定。
- **L5162-L15164（Elements MPMUpdatedLagrangian2D3N）**：将材料点分配到三节点背景单元的“挂靠”连接，用于形函数插值。质量与动量投影到背景节点的典型公式：
  - `[ m_I = \sum_p N_I(\mathbf{x}_p)\,m_p ]`
  - `[ \mathbf{p}_I = \sum_p N_I(\mathbf{x}_p)\,m_p\,\mathbf{v}_p ]`
- **L15165-20318（SubModelPart Parts_Parts_Auto1）**：收集属于材料域的节点/单元，便于求解器按名称加载。

## 背景网格：`source/granular_flow_2D_Grid.mdpa`
- **L1-L6**：`ModelPartData` 与空 `Properties`。
- **L7-L26136（Nodes）**：背景计算网格节点列表，定义欧拉网格用于求解动量方程。
- **L26137-67389（Elements Element2D3N）**：三角形背景单元，提供形函数 `N_I` 与梯度 `\nabla N_I`。网格动力平衡方程（弱式）离散为：
  - `[ \mathbf{M}\,\ddot{\mathbf{u}} + \mathbf{K}\,\mathbf{u} = \mathbf{f}_{\text{ext}} ]`
  其中质量矩阵 `[ \mathbf{M}_{IJ} = \int_{\Omega} \rho\,N_I N_J\,\mathrm{d}\Omega ]`，由上一步的粒子映射得到。
- **L67390-129982（SubModelParts）**：`Parts_Parts_Auto2` 包含网格节点与单元；`DISPLACEMENT_Displacement_Auto1` 子集用于施加零位移边界，对应前述约束过程。

## 运行流程速览（与 `simulation.Run()` 一致）
1. 解析 `ProjectParameters.json`，创建 `ModelPart` 并导入粒子/网格 (`L11-L38`)。
2. 读取 `ParticleMaterials.json` 绑定材料性质与本构 (`L1-L22`)。
3. 构建过程：施加重力、位移约束，准备输出。
4. 进入时间循环：每步将粒子量投影到背景网格，求解动量平衡获取 `[ \mathbf{a}_{n+1} ]`，用上文 Newmark 公式更新 `[ \mathbf{u}_{n+1}, \mathbf{v}_{n+1} ]`，再回写到材料点。
5. 按 `0.01 s` 输出 GiD/VTK 结果。

## 如何下载本说明
文档路径：`mpm/validation/granular_flow_2D/granular_flow_2D_explanation.md`。可直接下载原始文件：
`https://raw.githubusercontent.com/Lsq67opps/Kratos.Examples/copilot/create-markdown-explanation/mpm/validation/granular_flow_2D/granular_flow_2D_explanation.md`
