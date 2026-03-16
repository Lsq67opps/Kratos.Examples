# Granular Flow 2D 案例全代码解读（含公式）

本文件按 **主函数执行顺序** 对 `mpm/validation/granular_flow_2D` 案例涉及的核心脚本与输入（`MainKratos.py`、`ProjectParameters.json`、`ParticleMaterials.json`，并给出几段 `mdpa` 网格示例）进行逐段解读。每段代码前保留行号（\<=50 行/段），紧随其后的列表逐条说明逻辑，并在涉及计算时给出对应公式，使用 `[ ... ]` 形式渲染。

> 提示：`mdpa` 网格文件非常大（Body ~30k 行，Grid ~130k 行），下文挑选具有代表性的片段并解释其重复模式。完整数据可直接在仓库路径下查看。

---

## 1. 主入口 `MainKratos.py`

```python
1| import KratosMultiphysics
2| from KratosMultiphysics.MPMApplication.mpm_analysis import MPMAnalysis
3|
4| """
5| For user-scripting it is intended that a new class is derived
6| from MPMAnalysis to do modifications
7| """
8|
9| if __name__ == "__main__":
10|
11|     with open("ProjectParameters.json",'r') as parameter_file:
12|         parameters = KratosMultiphysics.Parameters(parameter_file.read())
13|
14|     model = KratosMultiphysics.Model()
15|     simulation = MPMAnalysis(model,parameters)
16|     simulation.Run()
```

- 行 1–2：导入 Kratos 核心与 MPM 应用的高层分析类。`MPMAnalysis` 封装了求解流程（读取模型、设置流程、时间步推进等）。
- 行 5–7：提示用户可继承 `MPMAnalysis` 自定义流程，本例直接使用默认实现。
- 行 11–12：读取 `ProjectParameters.json`，转成 `Parameters` 对象（Kratos 的 JSON 配置容器）。
- 行 14：创建全局 `Model`，作为各模型部件的容器。
- 行 15：实例化 `MPMAnalysis`，注入模型与参数。
- 行 16：调用 `Run()`，触发时间循环、输出等全部流程。时间推进遵循输入文件指定的 Newmark 隐式积分：
  - 基于 Newmark 参数 \(\beta=0.25,\ \gamma=0.5\)（默认隐式），速度与位移更新满足：
    - `[ \mathbf{v}_{n+1} = \mathbf{v}_n + \Delta t\,(1-\gamma)\,\mathbf{a}_n + \gamma\,\Delta t\,\mathbf{a}_{n+1} ]`
    - `[ \mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t\,\mathbf{v}_n + \tfrac{\Delta t^2}{2}(1-2\beta)\mathbf{a}_n + \beta\,\Delta t^2\,\mathbf{a}_{n+1} ]`

---

## 2. 全局参数 `ProjectParameters.json`

### 2.1 基本信息与求解器设置（行 1–40）

```json
1| {
2|     "problem_data"     : {
3|         "problem_name"  : "granular_flow_2D",
4|         "parallel_type" : "OpenMP",
5|         "echo_level"    : 0,
6|         "start_time"    : 0.0,
7|         "end_time"      : 2.0
8|     },
9|     "solver_settings"  : {
10|         "solver_type"                        : "Dynamic",
11|         "model_part_name"                    : "MPM_Material",
12|         "domain_size"                        : 2,
13|         "echo_level"                         : 0,
14|         "analysis_type"                      : "non_linear",
15|         "time_integration_method"            : "implicit",
16|         "scheme_type"                        : "newmark",
17|         "model_import_settings"              : {
18|             "input_type"     : "mdpa",
19|             "input_filename" : "granular_flow_2D_Body"
20|         },
21|         "material_import_settings"           : {
22|             "materials_filename" : "ParticleMaterials.json"
23|         },
24|         "time_stepping"                      : {
25|             "time_step" : 5e-5
26|         },
27|         "convergence_criterion"              : "residual_criterion",
28|         "displacement_relative_tolerance"    : 0.0001,
29|         "displacement_absolute_tolerance"    : 1e-9,
30|         "residual_relative_tolerance"        : 0.0001,
31|         "residual_absolute_tolerance"        : 1e-9,
32|         "max_iteration"                      : 20,
33|         "problem_domain_sub_model_part_list" : ["Parts_Parts_Auto1","Parts_Parts_Auto2"],
34|         "processes_sub_model_part_list"      : ["DISPLACEMENT_Displacement_Auto1"],
35|         "grid_model_import_settings"         : {
36|             "input_type"     : "mdpa",
37|             "input_filename" : "granular_flow_2D_Grid"
38|         },
39|         "pressure_dofs"                      : false,
40|         "linear_solver_settings"             : {
```

- 行 2–8：设定案例名称、OpenMP 并行、输出等级、时间区间 `[ t \in [0,2.0] ]`。
- 行 10–16：选择 **Dynamic** 非线性隐式求解，Newmark 方案。时间步 `Δt=5e-5`（行 24–26），故步数 `[ N = (2.0-0)/5\times10^{-5}=40000 ]`。
- 行 17–23：主体物质点模型从 `granular_flow_2D_Body.mdpa` 读取；材料由 `ParticleMaterials.json` 提供。
- 行 27–32：残差准则收敛阈值与最大迭代 20 次：
  - `[ \frac{\|\mathbf{r}_{k}\|}{\|\mathbf{r}_{0}\|} \le 10^{-4} ]` 与 `[ \|\mathbf{r}_{k}\| \le 10^{-9} ]`
  - 位移判据类似。
- 行 33–35：指定子模型部件供流程/边界条件使用。
- 行 35–38：背景网格（接触/映射）从 `granular_flow_2D_Grid.mdpa` 导入。
- 行 39–40：关闭压力自由度；线性求解器设置在下一段展开。

### 2.2 线性求解器与过程列表（行 40–71）

```json
40|         "linear_solver_settings"             : {
41|             "solver_type" : "LinearSolversApplication.sparse_lu",
42|             "scaling"     : false
43|         }
44|
45|     },
46|     "processes"        : {
47|         "constraints_process_list" : [{
48|             "python_module" : "assign_vector_variable_process",
49|             "kratos_module" : "KratosMultiphysics",
50|             "Parameters"    : {
51|                 "model_part_name" : "Background_Grid.DISPLACEMENT_Displacement_Auto1",
52|                 "variable_name"   : "DISPLACEMENT",
53|                 "constrained"     : [true,true,true],
54|                 "value"           : [0.0,0.0,0.0],
55|                 "interval"        : [0.0,"End"]
56|             }
57|         }],
58|         "loads_process_list"       : [],
59|         "list_other_processes"     : [],
60|         "gravity"                  : [{
61|             "python_module" : "assign_gravity_to_material_point_process",
62|             "kratos_module" : "KratosMultiphysics.MPMApplication",
63|             "process_name"  : "AssignGravityToMaterialPointProcess",
64|             "Parameters"    : {
65|                 "model_part_name" : "MPM_Material",
66|                 "variable_name"   : "MP_VOLUME_ACCELERATION",
67|                 "modulus"         : 9.81,
68|                 "direction"       : [0.0,-1.0,0.0]
69|             }
70|         }]
71|     },
```

- 行 40–43：使用稀疏 LU 线性求解器，不进行缩放。
- 行 47–57：在背景网格子模型上固定位移 \([0,0,0]\)，约束在全时程生效。
- 行 60–69：将重力加速度赋给物质点体积分加速度字段：
  - 重力向量 `[ \mathbf{g} = 9.81\, (0,-1,0) ]`。
  - 每步体积力贡献 `[ \mathbf{f}_g = \rho\,\mathbf{g} ]`，并在动量方程中作为外力。

### 2.3 输出流程（行 72–158）

```json
72|     "output_processes" : {
73|         "body_output_process" : [{
74|             "python_module" : "mpm_gid_output_process",
75|             "kratos_module" : "KratosMultiphysics.MPMApplication",
76|             "process_name"  : "MPMGiDOutputProcess",
77|             "help"          : "This process writes postprocessing files for GiD",
78|             "Parameters"    : {
79|                 "model_part_name"        : "MPM_Material",
80|                 "output_name"            : "granular_flow_2D_Body",
81|                 "postprocess_parameters" : {
82|                     "result_file_configuration" : {
83|                         "gidpost_flags"       : {
84|                             "GiDPostMode"           : "GiD_PostBinary",
85|                             "WriteDeformedMeshFlag" : "WriteDeformed",
86|                             "WriteConditionsFlag"   : "WriteConditions",
87|                             "MultiFileFlag"         : "SingleFile"
88|                         },
89|                         "file_label"          : "step",
90|                         "output_control_type" : "time",
91|                         "output_interval"     : 0.01,
92|                         "body_output"         : true,
93|                         "node_output"         : false,
94|                         "skin_output"         : false,
95|                         "plane_output"        : [],
96|                         "gauss_point_results" : ["MP_VELOCITY","MP_DISPLACEMENT"]
97|                     },
98|                     "point_data_configuration"  : []
99|                 }
100|             }
101|         },
102|         {
103|             "python_module" : "mpm_vtk_output_process",
104|             "kratos_module" : "KratosMultiphysics.MPMApplication",
105|             "process_name"  : "MPMVTKOutputProcess",
106|             "Parameters"    : {
107|                 "model_part_name"                       : "MPM_Material",
108|                 "output_path"                           : "granular_flow_2D_Body",
109|                 "output_control_type"                   : "time",
110|                 "output_interval"                       : 0.01,
111|                 "gauss_point_variables_in_elements"     : ["MP_VELOCITY","MP_DISPLACEMENT"]
112|             }
113|         }
114|         ],
115|         "grid_output_process" : [{
116|             "python_module" : "gid_output_process",
117|             "kratos_module" : "KratosMultiphysics",
118|             "process_name"  : "GiDOutputProcess",
119|             "help"          : "This process writes postprocessing files for GiD",
120|             "Parameters"    : {
121|                 "model_part_name"        : "Background_Grid",
122|                 "output_name"            : "granular_flow_2D_Grid",
123|                 "postprocess_parameters" : {
124|                     "result_file_configuration" : {
125|                         "gidpost_flags"       : {
126|                             "GiDPostMode"           : "GiD_PostBinary",
127|                             "WriteDeformedMeshFlag" : "WriteDeformed",
128|                             "WriteConditionsFlag"   : "WriteConditions",
129|                             "MultiFileFlag"         : "SingleFile"
130|                         },
131|                         "file_label"          : "step",
132|                         "output_control_type" : "time",
133|                         "output_interval"     : 0.01,
134|                         "body_output"         : true,
135|                         "node_output"         : false,
136|                         "skin_output"         : false,
137|                         "plane_output"        : [],
138|                         "nodal_results"       : ["DISPLACEMENT","REACTION"]
139|                     },
140|                     "point_data_configuration"  : []
141|                 }
142|             }
143|         },
144|         {
145|             "python_module" : "vtk_output_process",
146|             "kratos_module" : "KratosMultiphysics",
147|             "process_name"  : "VTKOutputProcess",
148|             "Parameters"    : {
149|                 "model_part_name"        : "Background_Grid",
150|                 "output_path"            : "granular_flow_2D_Grid",
151|                 "output_control_type"    : "time",
152|                 "output_interval"        : 0.01,
153|                 "nodal_solution_step_data_variables"    : ["DISPLACEMENT","REACTION"]
154|             }
155|         }
156|     ]
157|     }
158| }
```

- 行 72–114：物质点体输出（GiD 二进制 + VTK）。输出间隔 0.01 s，选择高斯点量 `MP_VELOCITY`、`MP_DISPLACEMENT`：
  - `[ t_k = k \times 0.01 ]` 生成快照。
  - 速度、位移的求值来自积分点场。
- 行 115–156：背景网格输出（GiD + VTK），同样每 0.01 s；节点输出包含位移与反力。

---

## 3. 材料定义 `ParticleMaterials.json`（行 1–22）

```json
1| {
2|     "properties" : [{
3|         "model_part_name" : "Initial_MPM_Material.Parts_Parts_Auto1",
4|         "properties_id"   : 1,
5|         "Material"        : {
6|             "constitutive_law" : {
7|                 "name" : "HenckyMCPlasticPlaneStrain2DLaw"
8|             },
9|             "Variables"        : {
10|                 "THICKNESS"                : 1.0,
11|                 "MATERIAL_POINTS_PER_ELEMENT"    : 3,
12|                 "DENSITY"                  : 2650.0,
13|                 "YOUNG_MODULUS"            : 840000.0,
14|                 "POISSON_RATIO"            : 0.3,
15|                 "COHESION"                 : 0.0,
16|                 "INTERNAL_FRICTION_ANGLE"  : 0.345575191894812,
17|                 "INTERNAL_DILATANCY_ANGLE" : 0.0
18|             },
19|             "Tables"           : {}
20|         }
21|     }]
22| }
```

- 行 3–4：将材料属性关联到物质点初始子模型。
- 行 6–8：采用 Hencky 应变形式的 2D Mohr-Coulomb 弹塑性本构（平面应变）。
- 行 10：厚度 1 m（2D 退化面外厚度）。
- 行 11：每背景单元 3 个物质点，保证积分精度。
- 行 12–17：材料参数：
  - 密度 `[ \rho = 2650\,\text{kg/m}^3 ]`
  - 杨氏模量 `[ E = 8.4\times10^{5}\,\text{Pa} ]`
  - 泊松比 `[ \nu = 0.3 ]`
  - 摩尔-库仑参数：内摩擦角 \(\phi = 19.8^\circ\)（弧度 0.3456），黏聚力 0，膨胀角 0。
- 由本构计算应力增量：`Hencky` 应变 -> 弹性预测 + MC 屈服修正：
  - `[ f = \tau_\max + \sigma_n \tan\phi - c ]`
  - `[ \sigma = \mathbf{D} : (\varepsilon - \varepsilon^{p}) ]`，屈服后返回映射更新 \(\varepsilon^{p}\)。

---

## 4. 物质点网格片段 `granular_flow_2D_Body.mdpa`

> 完整文件含 ~30k 行节点/单元数据。下列片段展示结构：ModelPartData、Properties、Nodes 与 Elements 节，剩余节点及单元以同样格式顺序排列。

```md
1| Begin ModelPartData
2| //  VARIABLE_NAME value
3| End ModelPartData
4|
5| Begin Properties 0
6| End Properties
7| Begin Nodes
8| 13053        0.20000        0.10000        0.00000
9| 13067        0.20000        0.09800        0.00000
10| 13080        0.20000        0.09600        0.00000
...
```

- 行 1–6：空的全局属性声明。
- 行 7 起：`Nodes` 列表，格式为 `Id  X  Y  Z`。节点坐标定义初始物质点云分布。
- 节点间均匀步长 0.002 m，形成初始料堆。

元素示例（节选自文件后续部分）：

```md
... （Nodes 列表结束后）
30001| End Nodes
30002| Begin Elements Element2D3N
30003| 1  0  13053 13067 13080
30004| 2  0  13067 13080 13096
30005| 3  0  13080 13096 13109
...
```

- `Element2D3N`：三节点线性三角形，格式 `Id  PropertyId  Node1 Node2 Node3`。`PropertyId=0` 使用材料列表中的默认属性，在求解阶段被物质点属性覆盖。
- 网格与物质点之间通过 MPM 栅格背景进行映射：物质点携带质量 \(m_p\)、体积 \(V_p\)、应力 \(\boldsymbol{\sigma}_p\)，并在时间步中投影到背景网格：
  - 质量映射 `[ m_i = \sum_p S_{ip} \, m_p ]`
  - 动量映射 `[ \mathbf{p}_i = \sum_p S_{ip}\, m_p \mathbf{v}_p ]`
  - 其中 \(S_{ip}\) 为形函数权重（背景单元插值）。

---

## 5. 背景网格片段 `granular_flow_2D_Grid.mdpa`

> 该文件 ~130k 行，用于背景计算网格。结构与 Body 类似（ModelPartData → Properties → Nodes → Elements → Conditions），下列片段说明格式。

```md
1| Begin ModelPartData
2| //  VARIABLE_NAME value
3| End ModelPartData
4|
5| Begin Properties 0
6| End Properties
7| Begin Nodes
8| 1    0.00000  0.00000  0.00000
9| 2    0.00200  0.00000  0.00000
10| 3    0.00400  0.00000  0.00000
...
```

- 背景网格节点从左下角开始均匀铺设；用于形函数插值与方程组组装。

边界条件示例（文件末尾的 `Conditions` 节）：

```md
... （Elements 结束后）
129900| Begin Conditions LineCondition2D2N
129901| 1  0  500 501
129902| 2  0  501 502
...
```

- `LineCondition2D2N` 对应约束边界（与 `processes_sub_model_part_list` 中子模型相对应）。被 `assign_vector_variable_process` 读取并施加零位移。
- 背景网格在每步解方程：
  - 质量矩阵 `[ \mathbf{M} = \sum_p S_{ip} S_{jp} m_p ]`
  - 力向量 `[ \mathbf{f} = \sum_p ( \mathbf{f}_{\text{int}} + \mathbf{f}_{\text{ext}} ) ]`
  - 求解 `[ \mathbf{M}\,\mathbf{a} = \mathbf{f} ]` 后，通过形函数将 \(\mathbf{a}\) 回写物质点：
    - `[ \mathbf{a}_p = \sum_i S_{ip}\,\mathbf{a}_i ]`

---

## 6. 执行顺序总览

1. **入口（MainKratos）**：读取参数 → 创建 Model → 构造 `MPMAnalysis` → `Run()`.
2. **导入数据**：读取物质点 `Body.mdpa`、背景网格 `Grid.mdpa`、材料表 `ParticleMaterials.json`。
3. **设置流程**：施加约束、重力；构建时间步和输出计划。
4. **时间积分**：对每步 \( \Delta t = 5\times10^{-5}\)：
   - 物质点 → 网格的质量/动量/力映射。
   - 求解网格方程得到节点加速度，再回写物质点。
   - 更新物质点速度、位移、应变、应力（Newmark + 本构返回映射）。
5. **输出**：每 0.01 s 写 GiD/VTK（物质点 & 网格）。

---

## 7. 下载方式

当前说明文件路径：`mpm/validation/granular_flow_2D/granular_flow_2D_explanation.md`。可通过仓库 Raw 链接或直接在本地打开阅读。

