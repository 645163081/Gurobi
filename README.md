Gurobi-EEDAPFSP Solver
This project provides a Gurobi-based Mixed Integer Linear Programming (MILP) implementation to solve an energy-efficient scheduling problem in Distributed Assembly Flowshops (EEDAPFSP). The proposed model integrates job sequencing, machine allocation, and three energy-saving strategies: shutdown, speed adjustment, and job postponement. The optimization objectives are to minimize the makespan and the total energy consumption (TEC).

🔧 File Overview
gantsolver.py: The main script that builds and solves the MILP model using the Gurobi solver. It also visualizes the resulting schedule as a Gantt chart.

I_8_5_2_2_1.txt: A sample input instance file.

I_20_5_2_4_7.txt: The industrial case used in Section 5.5 of the paper, representing a real-world engineering scenario.

📌 How to Use
Prepare your input instance in .txt format (e.g., I_8_5_2_2_1.txt) following the required structure.

Set the instance_file variable in gantsolver.py to point to your own .txt file path.

Run the script.
python gantsolver.py

📂 Output
After running the solver, a folder named results/ will be automatically generated. It contains:

Gantt chart (gantt.png): Visual representation of the optimized schedule, including job-machine assignments and assembly operations.

Scheduling results file: Includes the job sequences, energy-efficient operation data, and timing information.

📚 Case Study
The file I_20_5_2_4_7.txt represents a real-world automotive engine block manufacturing scenario.
