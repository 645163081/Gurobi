import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import time
import gurobipy as gp
from gurobipy import GRB

# 设置英文+负号兼容
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['axes.unicode_minus'] = False

def read_txt_data(file_path):
    """读取文本数据并解析为数据块"""
    data = []
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            current_block = []
            in_data_section = False
            for line in file:
                clean_line = line.strip()

                if not clean_line:
                    if current_block:
                        data.append(current_block)
                        current_block = []
                    in_data_section = False
                    continue

                if clean_line.isalpha() and not in_data_section:
                    if current_block:
                        data.append(current_block)
                    data.append([clean_line])
                    current_block = []
                    continue

                numbers = re.findall(r'-?\d+', clean_line)
                if numbers:
                    in_data_section = True
                    current_block.extend(map(int, numbers))

            if current_block:
                data.append(current_block)

        return data
    except Exception as e:
        print(f"读取文件失败: {e}")
        return None


def parse_instance_data(txt_data):
    """解析实例数据为结构化格式"""
    if not txt_data or len(txt_data) < 3:
        print("数据不足")
        return None

    # 基本参数 (n, m, R)
    if len(txt_data[0]) < 3:
        print("基本参数数据不完整")
        return None
    n = txt_data[0][0]  # 作业数量
    m = txt_data[0][1]  # 机器数量
    R = txt_data[0][2]  # 工厂数量
    print(f"基本参数: n={n}, m={m}, R={R}")

    # 处理时间数据
    pt = [[0] * n for _ in range(m)]
    block_idx = 1

    for i in range(n):
        if block_idx >= len(txt_data):
            print(f"错误：缺少作业{i + 1}的处理时间数据")
            return None

        block = txt_data[block_idx]
        if len(block) < 2 * m:
            block.extend([0] * (2 * m - len(block)))

        for j in range(m):
            pt[j][i] = block[2 * j + 1]

        block_idx += 1

    # 装配集数量
    H = 0
    while block_idx < len(txt_data):
        block = txt_data[block_idx]
        if isinstance(block, list) and block and isinstance(block[0], int):
            H = block[0]
            block_idx += 1
            break
        block_idx += 1

    if H == 0:
        print("警告：未找到装配集数量，使用默认值0")

    print(f"装配集数量 H={H}")

    # 装配时间数据
    p_A = [0] * (H + 1)
    for h in range(1, H + 1):
        if block_idx >= len(txt_data):
            print(f"警告：缺少装配集{h}的数据")
            break

        block = txt_data[block_idx]
        if len(block) < 2:
            ass_time = 0
        else:
            ass_id, ass_time = block[0], block[1]
            if ass_id != h:
                print(f"警告：装配集ID不匹配（期望{h}，实际{ass_id}）")

        p_A[h] = ass_time
        block_idx += 1

    # 作业归属信息
    prod_jobs = [[] for _ in range(H + 1)]
    for i in range(1, n + 1):
        if block_idx >= len(txt_data):
            prod_id = 1 if H > 0 else 0
            if 1 <= prod_id <= H:
                prod_jobs[prod_id].append(i)
            continue

        block = txt_data[block_idx]
        if len(block) < 2:
            job_id = i
            prod_id = 1 if H > 0 else 0
        else:
            job_id, prod_id = block[0], block[1]

        if job_id != i:
            print(f"警告：作业ID不匹配（期望{i}，实际{job_id}）")

        if 1 <= prod_id <= H:
            prod_jobs[prod_id].append(job_id)

        block_idx += 1

    return {
        'n': n,
        'm': m,
        'R': R,
        'H': H,
        'pt': pt,
        'p_A': p_A,
        'prod_jobs': prod_jobs
    }


def plot_gantt(result, instance_data, output_dir="results", filename="gantt.png"):
    """Draw Gantt chart for distributed assembly flowshop scheduling"""

    # Extract instance data
    rr = instance_data['R']  # number of factories
    m = instance_data['m']   # number of machines per factory
    t = instance_data['H']   # number of products
    prod_jobs = instance_data['prod_jobs']  # product-to-jobs mapping

    # Create job-to-product map
    job_to_product = {}
    all_jobs = set()
    for h in range(1, t + 1):
        for job in prod_jobs[h]:
            job_to_product[job] = h
            all_jobs.add(job)
    all_jobs = sorted(list(all_jobs))

    # Assign color per job
    color_map = plt.cm.get_cmap("tab20", len(all_jobs))
    job_colors = {job: color_map(i) for i, job in enumerate(all_jobs)}

    # Extract job sequences
    job_sequence = result['job_sequence']
    succ = {}
    for (k, j) in job_sequence:
        if k not in succ:
            succ[k] = []
        succ[k].append(j)

    # Build factory-specific job sequences
    factories_jobs = []
    start_jobs = succ.get(0, [])
    for start_job in start_jobs:
        factory_seq = []
        current = start_job
        while current != 0 and current in succ:
            factory_seq.append(current)
            current = succ[current][0]
        factories_jobs.append(factory_seq)

    # Create plot
    fig, ax = plt.subplots(figsize=(16, 8))

    # Set y-ticks
    y_ticks = []
    y_labels = []
    for f in range(rr):
        for i in range(m):
            y_pos = f * m + i
            y_ticks.append(y_pos)
            y_labels.append(f"F{f+1}-M{i+1}")
    y_ticks.append(rr * m)
    y_labels.append("Assembly")
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_ylim(-1, rr * m + 1)

    # Plot jobs per factory & machine
    for f, jobs_in_factory in enumerate(factories_jobs):
        for machine_i in range(m):
            for job in jobs_in_factory:
                start = result['all_variables'].get(f'Sij[{machine_i},{job}]', 0)
                end = result['all_variables'].get(f'Cij[{machine_i},{job}]', 0)
                if end - start > 0:
                    color = job_colors.get(job, 'gray')
                    rect = patches.Rectangle(
                        (start, f * m + machine_i - 0.4),
                        end - start,
                        0.8,
                        facecolor=color,
                        edgecolor='black',
                        alpha=0.8
                    )
                    ax.add_patch(rect)
                    ax.text(
                        (start + end) / 2, f * m + machine_i,
                        f'J{job}',
                        ha='center', va='center', color='black', fontsize=8
                    )

    # Plot assembly line
    prod_sequence = result['prod_sequence']
    succ_asm = {}
    for (l, s) in prod_sequence:
        if l not in succ_asm:
            succ_asm[l] = []
        succ_asm[l].append(s)

    asm_start = succ_asm[0][0] if 0 in succ_asm else 0
    asm_jobs = []
    current_asm = asm_start
    while current_asm != 0:
        asm_jobs.append(current_asm)
        current_asm = succ_asm.get(current_asm, [0])[0]

    for s in asm_jobs:
        end = result['all_variables'].get(f'CAs[{s}]', 0)
        start = end - instance_data['p_A'][s]
        duration = end - start
        job_ids = prod_jobs[s]
        first_job = job_ids[0] if job_ids else s
        color = job_colors.get(first_job, 'gray')

        rect = patches.Rectangle(
            (start, rr * m - 0.4),
            duration,
            0.8,
            facecolor=color,
            edgecolor='black',
            alpha=0.9
        )
        ax.add_patch(rect)
        ax.text(
            (start + end) / 2, rr * m,
            f'P{s}',
            ha='center', va='center', color='black', fontsize=8
        )

    # Axes and legend
    ax.set_xlim(0, result['makespan'] * 1.1)
    ax.set_xlabel("Time")
    ax.set_title("Distributed Assembly Flowshop Gantt Chart")

    # Legend (by job)
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color, label=f'J{job}')
        for job, color in job_colors.items()
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout()

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        plt.close()
    else:
        plt.show()

def solve_eedapfsp(instance_data, output_dir="results"):
    energy_history = []
    """求解节能调度问题"""
    # 提取实例数据
    n = instance_data['n']  # 作业数量
    m = instance_data['m']  # 机器数量
    rr = instance_data['R']  # 工厂数量（分布式工厂）
    t = instance_data['H']  # 产品数量（装配集数量）
    pt = instance_data['pt']  # 处理时间矩阵 (m x n)
    p_A = instance_data['p_A']  # 装配时间列表 (长度 t+1)
    prod_jobs = instance_data['prod_jobs']  # 作业归属信息

    # 速度调整相关参数
    speeds = instance_data.get('speeds', [0.76, 1.0, 1.3])
    D = len(speeds)
    base_speed = instance_data.get('base_speed', 1.0)
    w1 = instance_data.get('w1', 0.5)

    # 空闲关停相关参数
    idle_power = instance_data.get('idle_power', 1.0)  # 空闲功率，单位kW
    shutdown_power = instance_data.get('shutdown_power', 4.5)  # 单次关停能耗
    shutdown_threshold = instance_data.get('shutdown_threshold', 4.50)  # 超过该空闲时间视为关停

    max_pt = max(max(row) for row in pt) / min(speeds)
    M = (n * m * max_pt) * 2

    # 创建模型
    model = gp.Model("EEDAPFSP")
    model.setParam('NumericFocus', 3)  # 提高数值精度

    # 添加变量
    X = model.addVars(n + 1, n + 1, vtype=GRB.BINARY, name="X")
    Y = model.addVars(t + 1, t + 1, vtype=GRB.BINARY, name="Y")
    Cij = model.addVars(m, n + 1, vtype=GRB.CONTINUOUS, name="Cij")
    CAs = model.addVars(t + 1, vtype=GRB.CONTINUOUS, name="CAs")
    Cmax = model.addVar(vtype=GRB.CONTINUOUS, name="Cmax")
    processing_energy = model.addVar(vtype=GRB.CONTINUOUS, name="processing_energy")
    AE = model.addVar(vtype=GRB.CONTINUOUS, name="AE")
    z = model.addVars(m, n + 1, D, vtype=GRB.BINARY, name="V")
    pd = model.addVars(m, n + 1, vtype=GRB.CONTINUOUS, name="pd")
    idle_time = model.addVars(
        [(i, j, k) for i in range(m)
         for j in range(1, n + 1) for k in range(1, n + 1)
         if j != k],
        vtype=GRB.CONTINUOUS, name="IdleTime"
    )
    is_shutdown = model.addVars(
        [(i, j, k) for i in range(m) for j in range(1, n + 1) for k in range(1, n + 1) if j != k],
        vtype=GRB.BINARY, name="IsShutdown"
    )
    idle_energy = model.addVar(vtype=GRB.CONTINUOUS, name="idle_energy")
    shutdown_energy = model.addVar(vtype=GRB.CONTINUOUS, name="shutdown_energy")
    Sij = model.addVars(m, n + 1, vtype=GRB.CONTINUOUS, name="Sij")  # 开始时间

    # ================== 基本约束 ================== #
    # 约束1-5：作业调度逻辑
    model.addConstrs((gp.quicksum(X[k, j] for k in range(n + 1) if k != j) == 1 for j in range(1, n + 1)), "constr1")
    model.addConstrs((gp.quicksum(X[k, j] for j in range(n + 1) if k != j) <= 1 for k in range(1, n + 1)), "constr2")
    model.addConstr(gp.quicksum(X[0, j] for j in range(1, n + 1)) == rr, "constr3")
    model.addConstr(gp.quicksum(X[k, 0] for k in range(1, n + 1)) == rr, "constr4")  # 修改：应为r个工厂的结束
    model.addConstrs((X[k, j] + X[j, k] <= 1 for j in range(1, n) for k in range(1, j)), "constr5")

    # 速度选择约束
    model.addConstrs((gp.quicksum(z[i, j, d] for d in range(D)) == 1 for i in range(m) for j in range(1, n + 1)))
    model.addConstrs((pd[i, j] == gp.quicksum(pt[i][j - 1] * (base_speed / speeds[d]) * z[i, j, d] for d in range(D))
                      for i in range(m) for j in range(1, n + 1)))
    model.addConstrs((pd[i, 0] == 0 for i in range(m)), "virtual_job_pd")

    # ================== 时间约束 ================== #
    model.addConstrs(
        (Sij[i, j] == Cij[i, j] - pd[i, j]
         for i in range(m) for j in range(n + 1)),
        name="start_time_def"
    )

    # 机器顺序约束
    for i in range(m):
        if i == 0:
            for j in range(1, n + 1):
                model.addConstr(Cij[i, j] >= pd[i, j], f"constr6_m{i}_j{j}")
        else:
            for j in range(1, n + 1):
                model.addConstr(Cij[i, j] >= Cij[i - 1, j] + pd[i, j], f"constr6_m{i}_j{j}")

    # 作业顺序约束 - 使用指标约束确保数值稳定性
    for i in range(m):
        for j in range(1, n + 1):
            for k in range(n + 1):
                if k != j:
                    model.addGenConstrIndicator(
                        X[k, j], True,
                        Cij[i, j] >= Cij[i, k] + pd[i, j],
                        name=f"constr7_m{i}_k{k}_j{j}"
                    )

    # 虚拟作业约束
    for i in range(m):
        for j in range(1, n + 1):
            model.addGenConstrIndicator(
                X[0, j], True,
                Cij[i, j] >= pd[i, j],
                name=f"constr7_virtual_m{i}_j{j}"
            )

    # ================== 装配约束 ================== #
    model.addConstrs((gp.quicksum(Y[l, s] for l in range(t + 1) if l != s) == 1 for s in range(1, t + 1)), "constr8")
    model.addConstrs((gp.quicksum(Y[l, s] for s in range(t + 1) if l != s) <= 1 for l in range(t + 1)), "constr9")
    model.addConstrs((Y[l, s] + Y[s, l] <= 1 for s in range(1, t + 1) for l in range(1, s)), "constr10")

    # 作业-产品关系
    G = [[0] * (t + 1) for _ in range(n + 1)]
    for h in range(1, t + 1):
        for job in prod_jobs[h]:
            if job <= n:
                G[job][h] = 1

    # 装配前完成约束
    for s in range(1, t + 1):
        for j in range(1, n + 1):
            if G[j][s] == 1:
                model.addConstr(
                    CAs[s] >= Cij[m - 1, j] + p_A[s],
                    f"constr11_s{s}_j{j}"
                )

    # 装配顺序约束
    for s in range(1, t + 1):
        for l in range(t + 1):
            if l != s:
                model.addGenConstrIndicator(
                    Y[l, s], True,
                    CAs[s] >= CAs[l] + p_A[s],
                    name=f"constr12_l{l}_s{s}"
                )

    # ================== 空闲及关停能耗约束 ================== #
    for i in range(m):
        for j in range(1, n + 1):
            for k in range(1, n + 1):
                if j != k and j != 0:
                    model.addGenConstrIndicator(
                        X[j, k], True,
                        idle_time[i, j, k] == Sij[i, k] - Cij[i, j],
                        name=f"idle_time_exact_m{i}_j{j}_k{k}"
                    )
                    model.addGenConstrIndicator(
                        X[j, k], False,
                        idle_time[i, j, k] == 0,
                        name=f"idle_zero_indicator_m{i}_j{j}_k{k}"
                    )
                    model.addConstr(idle_time[i, j, k] >= 0, f"idle_time_nonneg_m{i}_j{j}_k{k}")



    # ================== 目标函数相关约束 ================== #
    model.addConstrs((Cmax >= CAs[s] for s in range(1, t + 1)), "constr13")
    model.addConstrs((Cij[i, 0] == 0 for i in range(m)), "virtual_job0")
    model.addConstr(CAs[0] == 0, "virtual_product0")

    # 加工能耗
    model.addConstr(
        processing_energy == gp.quicksum(
            gp.quicksum(
                gp.quicksum(
                    4 * speeds[d] * pt[i][j - 1] * z[i, j, d]
                    for d in range(D)
                )
                for j in range(1, n + 1)
            )
            for i in range(m)
        ),
        name="processing_energy_def"
    )

    # 关停能耗
    model.addConstr(
        shutdown_energy == gp.quicksum(
            shutdown_power * is_shutdown[i, j, k]
            for i in range(m) for j in range(1, n + 1) for k in range(1, n + 1) if j != k
        ),
        name="shutdown_energy_def"
    )

    # 空闲能耗
    model.addConstr(
        idle_energy == gp.quicksum(
            idle_time[i, j, k] * idle_power * (1 - is_shutdown[i, j, k])
            for i in range(m) for j in range(1, n + 1) for k in range(1, n + 1)
            if j != k
        ),
        name="idle_energy_def"
    )

    model.addConstr(AE == processing_energy + idle_energy + shutdown_energy, name="AE_total")

    # ================== 目标函数 ================== #
    model.setObjective(
        w1 * Cmax + (1 - w1) * AE,
        GRB.MINIMIZE
    )
    model.setParam('TimeLimit', 600)

    # 求解模型
    start_time = time.time()
    # 添加能耗记录列表
    energy_history = []
    # 回调函数记录目标值变化
    def energy_callback(model, where):
        if where == GRB.Callback.MIP:
            try:
                obj_val = model.cbGet(GRB.Callback.MIP_OBJBST)
                if obj_val < GRB.INFINITY:
                    energy_history.append(obj_val)
            except Exception:
                pass

    # 启动优化并记录能耗变化
    model.optimize(energy_callback)
    solve_time = time.time() - start_time

    # 准备结果
    result = {
        'status': model.status,
        'solve_time': solve_time
    }

    if model.status == GRB.OPTIMAL:
        result['makespan'] = Cmax.X
        result['energy'] = AE.X
        result['weighted_objective'] = model.ObjVal
        result['processing_energy'] = processing_energy.X
        result['idle_energy'] = idle_energy.X
        result['shutdown_energy'] = shutdown_energy.X

        # 输出所有变量值
        all_vars = {}
        for v in model.getVars():
            all_vars[v.VarName] = v.X
        result['all_variables'] = all_vars

        # 提取调度信息
        job_sequence = [(k, j) for k in range(n + 1) for j in range(n + 1) if k != j and X[k, j].X > 0.5]
        prod_sequence = [(l, s) for l in range(t + 1) for s in range(t + 1) if l != s and Y[l, s].X > 0.5]
        speed_selection = {(i, j): (d, speeds[d], pd[i, j].X)
                           for i in range(m) for j in range(1, n + 1) for d in range(D) if z[i, j, d].X > 0.5}

        result['job_sequence'] = job_sequence
        result['prod_sequence'] = prod_sequence
        result['speed_selection'] = speed_selection

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            with open(os.path.join(output_dir, "result.txt"), 'w') as f:
                f.write(f"最优makespan: {Cmax.X}\n总能耗: {AE.X}\n\n能耗分项：\n")
                f.write(
                    f"加工能耗: {processing_energy.X:.2f}\n空闲能耗: {idle_energy.X:.2f}\n关停能耗: {shutdown_energy.X:.2f}\n")
                f.write(f"加权目标值: {model.ObjVal}\n求解时间: {solve_time:.2f}秒\n\n权重设置: w1 = {w1}\n\n")

                f.write("作业调度序列:\n")
                for k, j in job_sequence:
                    f.write(f"作业{k} → 作业{j}\n")

                f.write("\n产品调度序列:\n")
                for l, s in prod_sequence:
                    f.write(f"产品{l} → 产品{s}\n")

                f.write("\n作业完成时间:\n")
                for j in range(n + 1):
                    f.write(f"作业{j}: ")
                    for i in range(m):
                        f.write(f"机器{i}={Cij[i, j].X:.1f} ")
                    f.write("\n")

                f.write("\n作业开始时间:\n")
                for j in range(n + 1):
                    f.write(f"作业{j}: ")
                    for i in range(m):
                        f.write(f"机器{i}={Sij[i, j].X:.1f} ")
                    f.write("\n")

                f.write("\n产品完成时间:\n")
                for s in range(t + 1):
                    f.write(f"产品{s}: {CAs[s].X:.1f}\n")

                f.write("\n速度选择:\n")
                for i in range(m):
                    f.write(f"机器 {i}:\n")
                    for j in range(1, n + 1):
                        if (i, j) in speed_selection:
                            d, speed, proc_time = speed_selection[(i, j)]
                            f.write(f"  作业 {j}: 速度等级 {d} (速度={speed}), 加工时间={proc_time:.1f}\n")
                        else:
                            f.write(f"  作业 {j}: 未选择速度 (虚拟作业)\n")

                f.write("\n延后分析:\n")
                for j in range(1, n + 1):
                    for i in range(m):
                        earliest_start = 0
                        # 计算最早可能开始时间
                        if i > 0:
                            earliest_start = max(earliest_start, Cij[i - 1, j].X)
                        # 考虑前序作业
                        for k in range(1, n + 1):
                            if X[k, j].X > 0.5:  # k是j的前序
                                earliest_start = max(earliest_start, Cij[i, k].X)

                        actual_start = Sij[i, j].X
                        delay = actual_start - earliest_start
                        if delay > 1e-5:
                            f.write(
                                f"作业{j}在机器{i}上延后了: {delay:.2f} (最早: {earliest_start:.2f}, 实际: {actual_start:.2f})\n")

        # 绘制甘特图
        try:
            plot_gantt(result, instance_data, output_dir)
        except Exception as e:
            print(f"绘制甘特图时出错: {e}")

    return result

def run_single_instance(file_path):
    """运行单个实例并保存结果"""
    # 创建输出目录
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'=' * 50}")
    print(f"处理实例: {file_path}")

    # 检查文件是否存在
    if not os.path.exists(file_path):
        print(f"文件不存在: {file_path}")
        return

    # 读取实例数据
    txt_data = read_txt_data(file_path)
    if txt_data is None:
        print(f"跳过实例: {file_path} (读取失败)")
        return

    # 解析实例数据
    instance_data = parse_instance_data(txt_data)
    if instance_data is None:
        print(f"跳过实例: {file_path} (解析失败)")
        return

    print(f"实例参数: 作业={instance_data['n']}, 机器={instance_data['m']}, "
          f"工厂={instance_data['R']}, 装配集={instance_data['H']}")

    # 提取实例名称
    file_name = os.path.basename(file_path)
    instance_name = os.path.splitext(file_name)[0]

    # 求解实例
    instance_dir = os.path.join(output_dir, instance_name)
    result = solve_eedapfsp(instance_data, instance_dir)

# 使用示例
if __name__ == "__main__":
    # 指定要处理的实例文件路径
    instance_file = "D:/project/研究生/Gurobi/I_8_5 _2_2_1.txt"

    run_single_instance(instance_file)