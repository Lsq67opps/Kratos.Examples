# granular_flow_2D 案例逐行解读（按主函数运行顺序）

> 说明：本文件按执行顺序对 `mpm/validation/granular_flow_2D/source` 下的脚本与配置逐行/分段解释。每段先复述代码（含行号），再给出对应解释与计算公式。公式均以 `[ ... ]` 形式渲染。若后续内容超过 1000 行，将在新文件继续，但当前内容不足 1000 行。

---

## 1. `source/MainKratos.py`（入口脚本）

### 1.1 代码（含行号）
```text
L1 : import KratosMultiphysics
L2 : from KratosMultiphysics.MPMApplication.mpm_analysis import MPMAnalysis

L5 : """
L6 : For user-scripting it is intended that a new class is derived
L7 : from MPMAnalysis to do modifications
L8 : """

L9 : if __name__ == "__main__":
L11:     with open("ProjectParameters.json",'r') as parameter_file:
L12:         parameters = KratosMultiphysics.Parameters(parameter_file.read())
L14:     model = KratosMultiphysics.Model()
L15:     simulation = MPMAnalysis(model,parameters)
L16:     simulation.Run()
```

### 1.2 解释与公式
- **L1-L2**：导入 Kratos 核心和 MPM 求解器封装。`MPMAnalysis` 提供初始化、时间步循环、输出等高层流程。
- **L5-L8**：文档字符串，提示可派生自 `MPMAnalysis` 定制。
- **L9-L16**：主流程：
  1. 读取 `ProjectParameters.json` 构造 `Parameters`。
  2. 创建空 `Model` 容器。
  3. 用参数实例化 `MPMAnalysis`。
  4. 调用 `simulation.Run()`，内部执行：导入模型、读材料、施加过程、时间积分、输出。
- 时间积分（Newmark 隐式，Kratos 默认 `beta, gamma`）：  
  `[ \mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t\,\mathbf{v}_n + \Delta t^2\!\left(\tfrac{1}{2}-\beta\right)\mathbf{a}_n + \beta\,\Delta t^2\,\mathbf{a}_{n+1} ]`  
  `[ \mathbf{v}_{n+1} = \mathbf{v}_n + \Delta t\!\left((1-\gamma)\mathbf{a}_n + \gamma\,\mathbf{a}_{n+1}\right) ]`

---

## 2. `source/ProjectParameters.json`（全局与求解设置）

### 2.1 代码（含行号）
```text
L2 : "problem_data" : {
L3 :     "problem_name"  : "granular_flow_2D",
L4 :     "parallel_type" : "OpenMP",
L5 :     "echo_level"    : 0,
L6 :     "start_time"    : 0.0,
L7 :     "end_time"      : 2.0
L8 : },
L9 : "solver_settings" : {
L10:     "solver_type"             : "Dynamic",
L11:     "model_part_name"         : "MPM_Material",
L12:     "domain_size"             : 2,
L13:     "echo_level"              : 0,
L14:     "analysis_type"           : "non_linear",
L15:     "time_integration_method" : "implicit",
L16:     "scheme_type"             : "newmark",
L17:     "model_import_settings"   : {
L18:         "input_type"     : "mdpa",
L19:         "input_filename" : "granular_flow_2D_Body"
L20:     },
L21:     "material_import_settings": {
L22:         "materials_filename" : "ParticleMaterials.json"
L23:     },
L24:     "time_stepping" : {
L25:         "time_step" : 5e-5
L26:     },
L27:     "convergence_criterion"           : "residual_criterion",
L28:     "displacement_relative_tolerance": 0.0001,
L29:     "displacement_absolute_tolerance": 1e-9,
L30:     "residual_relative_tolerance"     : 0.0001,
L31:     "residual_absolute_tolerance"     : 1e-9,
L32:     "max_iteration"                   : 20,
L33:     "problem_domain_sub_model_part_list" : ["Parts_Parts_Auto1","Parts_Parts_Auto2"],
L34:     "processes_sub_model_part_list"      : ["DISPLACEMENT_Displacement_Auto1"],
L35:     "grid_model_import_settings" : {
L36:         "input_type"     : "mdpa",
L37:         "input_filename" : "granular_flow_2D_Grid"
L38:     },
L39:     "pressure_dofs" : false,
L40:     "linear_solver_settings" : {
L41:         "solver_type" : "LinearSolversApplication.sparse_lu",
L42:         "scaling"     : false
L43:     }
L44: },
L47: "processes" : {
L48:     "constraints_process_list" : [{
L49:         "python_module" : "assign_vector_variable_process",
L50:         "kratos_module" : "KratosMultiphysics",
L51:         "Parameters"    : {
L52:             "model_part_name" : "Background_Grid.DISPLACEMENT_Displacement_Auto1",
L53:             "variable_name"   : "DISPLACEMENT",
L54:             "constrained"     : [true,true,true],
L55:             "value"           : [0.0,0.0,0.0],
L56:             "interval"        : [0.0,"End"]
L57:         }
L58:     }],
L59:     "loads_process_list"       : [],
L60:     "list_other_processes"     : [],
L61:     "gravity" : [{
L62:         "python_module" : "assign_gravity_to_material_point_process",
L63:         "kratos_module" : "KratosMultiphysics.MPMApplication",
L64:         "process_name"  : "AssignGravityToMaterialPointProcess",
L65:         "Parameters"    : {
L66:             "model_part_name" : "MPM_Material",
L67:             "variable_name"   : "MP_VOLUME_ACCELERATION",
L68:             "modulus"         : 9.81,
L69:             "direction"       : [0.0,-1.0,0.0]
L70:         }
L71:     }]
L72: },
L73: "output_processes" : {
L74:     "body_output_process" : [{
L75:         "python_module" : "mpm_gid_output_process",
L76:         "kratos_module" : "KratosMultiphysics.MPMApplication",
L77:         "process_name"  : "MPMGiDOutputProcess",
L79:         "Parameters"    : {
L80:             "model_part_name"        : "MPM_Material",
L81:             "output_name"            : "granular_flow_2D_Body",
L82:             "postprocess_parameters" : {
L83:                 "result_file_configuration" : {
L84:                     "gidpost_flags" : {
L85:                         "GiDPostMode"           : "GiD_PostBinary",
L86:                         "WriteDeformedMeshFlag" : "WriteDeformed",
L87:                         "WriteConditionsFlag"   : "WriteConditions",
L88:                         "MultiFileFlag"         : "SingleFile"
L89:                     },
L90:                     "file_label"          : "step",
L91:                     "output_control_type" : "time",
L92:                     "output_interval"     : 0.01,
L93:                     "body_output"         : true,
L94:                     "node_output"         : false,
L95:                     "skin_output"         : false,
L96:                     "plane_output"        : [],
L97:                     "gauss_point_results" : ["MP_VELOCITY","MP_DISPLACEMENT"]
L98:                 },
L99:                 "point_data_configuration" : []
L100:            }
L101:        }
L102:    },
L103:    {
L104:        "python_module" : "mpm_vtk_output_process",
L105:        "kratos_module" : "KratosMultiphysics.MPMApplication",
L106:        "process_name"  : "MPMVTKOutputProcess",
L107:        "Parameters"    : {
L108:            "model_part_name"                   : "MPM_Material",
L109:            "output_path"                       : "granular_flow_2D_Body",
L110:            "output_control_type"               : "time",
L111:            "output_interval"                   : 0.01,
L112:            "gauss_point_variables_in_elements" : ["MP_VELOCITY","MP_DISPLACEMENT"]
L113:        }
L114:    }],
L115:    "grid_output_process" : [{
L116:        "python_module" : "gid_output_process",
L117:        "kratos_module" : "KratosMultiphysics",
L118:        "process_name"  : "GiDOutputProcess",
L120:        "Parameters"    : {
L121:            "model_part_name"        : "Background_Grid",
L122:            "output_name"            : "granular_flow_2D_Grid",
L123:            "postprocess_parameters" : {
L124:                "result_file_configuration" : {
L125:                    "gidpost_flags" : {
L126:                        "GiDPostMode"           : "GiD_PostBinary",
L127:                        "WriteDeformedMeshFlag" : "WriteDeformed",
L128:                        "WriteConditionsFlag"   : "WriteConditions",
L129:                        "MultiFileFlag"         : "SingleFile"
L130:                    },
L131:                    "file_label"          : "step",
L132:                    "output_control_type" : "time",
L133:                    "output_interval"     : 0.01,
L134:                    "body_output"         : true,
L135:                    "node_output"         : false,
L136:                    "skin_output"         : false,
L137:                    "plane_output"        : [],
L138:                    "nodal_results"       : ["DISPLACEMENT","REACTION"]
L139:                },
L140:                "point_data_configuration" : []
L141:            }
L142:        }
L143:    },
L144:    {
L145:        "python_module" : "vtk_output_process",
L146:        "kratos_module" : "KratosMultiphysics",
L147:        "process_name"  : "VTKOutputProcess",
L148:        "Parameters"    : {
L149:            "model_part_name"     : "Background_Grid",
L150:            "output_path"         : "granular_flow_2D_Grid",
L151:            "output_control_type" : "time",
L152:            "output_interval"     : 0.01,
L153:            "nodal_solution_step_data_variables" : ["DISPLACEMENT","REACTION"]
L154:        }
L155:    }]
L156: }
```

### 2.2 解释与公式
- **L2-L8（problem_data）**：问题名、并行模式、时间区间 `[0,2]` 秒。时间步 `Δt = 5×10^{-5}`（见 L24-L26），步数 `[ N_{\text{step}} = \frac{2.0-0.0}{5\times10^{-5}} ]`。
- **L10-L16**：动态非线性、隐式 Newmark。`domain_size=2` 指 2D。
- **L17-L23 / L35-L38**：导入材料点与背景网格 mdpa 文件。
- **L24-L26（time_stepping）**：时间步长 `Δt`；进入 Newmark 公式同 §1.2。
- **L27-L32**：残差收敛准则与容差，迭代上限 20。
- **L33-L34**：子模型部件列表，供求解器和过程引用。
- **L39-L43**：稀疏 LU 线性解法。
- **L48-L57（约束过程）**：背景网格边界子模型 `DISPLACEMENT_Displacement_Auto1` 全约束：  
  `[ \mathbf{u} = (0,0,0) ]` 在整个时间区间。
- **L61-L70（重力）**：对材料点施加体加速度：  
  `[ \mathbf{a}_g = 9.81\,(0,-1,0) ]`。  
  动量平衡（连续体弱式简化）：`[ \rho\,\mathbf{a} = \nabla\!\cdot\!\boldsymbol{\sigma} + \rho\,\mathbf{a}_g ]`。
- **L74-L114（材料点输出）**：GiD 与 VTK 按时间步 0.01 s 输出 `MP_VELOCITY, MP_DISPLACEMENT`。
- **L115-L155（网格输出）**：GiD/VTK 输出背景网格的 `DISPLACEMENT, REACTION`，同样 0.01 s 间隔。

---

## 3. `source/ParticleMaterials.json`（材料与本构）

### 3.1 代码（含行号）
```text
L1 : {
L2 :     "properties" : [{
L3 :         "model_part_name" : "Initial_MPM_Material.Parts_Parts_Auto1",
L4 :         "properties_id"   : 1,
L5 :         "Material"        : {
L6 :             "constitutive_law" : {
L7 :                 "name" : "HenckyMCPlasticPlaneStrain2DLaw"
L8 :             },
L9 :             "Variables" : {
L10:                 "THICKNESS"                : 1.0,
L11:                 "MATERIAL_POINTS_PER_ELEMENT" : 3,
L12:                 "DENSITY"                  : 2650.0,
L13:                 "YOUNG_MODULUS"            : 840000.0,
L14:                 "POISSON_RATIO"            : 0.3,
L15:                 "COHESION"                 : 0.0,
L16:                 "INTERNAL_FRICTION_ANGLE"  : 0.345575191894812,
L17:                 "INTERNAL_DILATANCY_ANGLE" : 0.0
L18:             },
L19:             "Tables" : {}
L20:         }
L21:     }]
L22: }
```

### 3.2 解释与公式
- **L3-L4**：属性集绑定到材料点子模型 `Parts_Parts_Auto1`。
- **L6-L8**：本构采用平面应变 Hencky-Mohr-Coulomb 塑性。
- **L10-L17（变量）**：
  - 厚度 `[ t = 1.0 ]`
  - 每单元材料点数 3
  - 密度 `[ \rho = 2650\,\text{kg/m}^3 ]`
  - 弹性：`[ \boldsymbol{\sigma} = \mathbf{C}:\boldsymbol{\varepsilon} ]`，其中 `C` 由 `[ E=840000,\ \nu=0.3 ]` 生成的平面应变刚度矩阵。
  - 屈服函数（Mohr-Coulomb，粘聚力 0）：  
    `[ F = \sqrt{3J_2} + p\,\tan\varphi - c = 0 ]`，`[ \varphi = 0.3456\ \text{rad},\ c=0 ]`。
- **L19**：无额外材料表。

---

## 4. 数据文件说明（不逐行展开）
- `granular_flow_2D_Body.mdpa`：含材料点节点、MPMUpdatedLagrangian2D3N 连接、子模型部件。质量与动量映射到背景网格的典型关系：  
  `[ m_I = \sum_p N_I(\mathbf{x}_p)\,m_p ]`，`[ \mathbf{p}_I = \sum_p N_I(\mathbf{x}_p)\,m_p\,\mathbf{v}_p ]`。
- `granular_flow_2D_Grid.mdpa`：背景网格节点与三角形单元 Element2D3N，子模型包含边界约束集。网格动量方程离散：  
  `[ \mathbf{M}\,\ddot{\mathbf{u}} + \mathbf{K}\,\mathbf{u} = \mathbf{f}_{\text{ext}} ]`，其中 `[ \mathbf{M}_{IJ} = \int \rho N_I N_J\,\mathrm{d}\Omega ]` 由粒子映射得到。

---

## 5. 执行流程（与 `simulation.Run()` 对应）
1. 解析 `ProjectParameters.json`（§2），创建模型并导入粒子/网格 mdpa。
2. 读取 `ParticleMaterials.json`（§3）绑定本构与属性。
3. 构建过程：施加重力、位移约束；准备 GiD/VTK 输出。
4. 时间步循环：粒子 → 网格映射质量/动量，求解动量平衡得 `[ \mathbf{a}_{n+1} ]`，用 Newmark（§1.2）更新 `[ \mathbf{u}_{n+1}, \mathbf{v}_{n+1} ]`，再回写粒子。
5. 每 0.01 s 输出结果。

---

## 6. 下载
可直接下载本说明（raw）：  
`https://raw.githubusercontent.com/Lsq67opps/Kratos.Examples/copilot/create-markdown-explanation/mpm/validation/granular_flow_2D/granular_flow_2D_explanation.md`
