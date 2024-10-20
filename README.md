# Accelerating Distributed Optimization via Over-the-Air Computing

## Overview
This repository presents a novel framework for accelerating distributed optimization using Over-the-Air Computing (AirComp). Distributed optimization is essential for various applications like smart grid management, wireless resource allocation, machine learning, and large-scale network control. A significant challenge in such systems is the high communication overhead, which can become a bottleneck. To address this, we propose a method that leverages AirComp to achieve low-latency, efficient data aggregation in large-scale networks.

## Key Features

### Over-the-Air Computing for Distributed Optimization:
- AirComp exploits the superposition property of wireless multiple access channels (MAC) to perform simultaneous data aggregation, reducing the need for individual transmissions, and significantly cutting down bandwidth usage.
  
### Distributed Primal-Dual Subgradient Optimization:
- The framework employs a Distributed Primal-Dual (DPD) subgradient method to handle general convex optimization problems. Each device computes updates locally and transmits them via AirComp, enabling fast global updates without extensive data exchanges.

### Robust to Channel Fading and Noise:
- Despite the presence of channel fading and additive noise, the proposed method ensures that the optimization constraints are met and the problem remains feasible over time. Proper power control helps mitigate the effects of noise and achieve near-optimal performance.

## System Model
The system consists of a wireless network where multiple devices work together with a central server to solve a global optimization problem. Each device handles a part of the problem locally, and their solutions are combined through AirComp to form a global solution.

Key components:
- **Devices**: Perform local optimization steps and transmit updates via AirComp.
- **Central Server**: Aggregates the transmitted updates to adjust the global variables and broadcasts them back to the devices for further refinement.
- **OTA**: Utilizes the superposition of wireless signals to achieve efficient data aggregation, leading to significant bandwidth savings compared to traditional methods.

### Two Practical Use Cases:
1. **Smart Grid Energy Management**:
   - The DPD-AirComp framework is applied to optimize energy distribution in a smart grid, ensuring efficient energy use and grid stability. A Stackelberg game is used to model the interaction between the smart energy manager and electric vehicles (EVs), where the manager sets prices and the EVs adjust their energy demands based on pricing.

2. **FDMA Wireless Resource Allocation**:
   - The framework is used to allocate power and bandwidth in a frequency-division multiple access (FDMA) system. The goal is to maximize the sum-rate across multiple users while maintaining quality-of-service (QoS) constraints. The DPD-AirComp framework optimizes resource allocation in real time, providing fast convergence compared to traditional methods.

## Theoretical Guarantees
- **Zero Constraint Violation**: The proposed DPD-AirComp algorithm achieves zero expected constraint violation as the number of iterations grows, ensuring feasibility of the optimization problem even in the presence of partial user participation.
  
- **Optimality Gap**: By using power control, the framework minimizes the expected optimality gap, ensuring that the solution closely approximates the true optimal solution.

## Results
- The proposed DPD-AirComp framework achieves **an order of magnitude faster convergence** compared to traditional digital orthogonal multiple access methods (e.g., TDMA).
- **Smart Grid Management**: Significant improvements in energy distribution efficiency are achieved, with fast convergence to optimal pricing and energy allocation.
- **Wireless Resource Allocation**: DPD-AirComp enhances resource utilization, providing near-optimal solutions with substantially lower communication overhead.

## Configuration
The system parameters can be adjusted depending on the application scenario:
- Number of devices
- Channel fading and noise levels
- Power constraints

## References
This work is based on the following paper:

N. A. Mitsiou, P. S. Bouzinis, P. D. Diamantoulakis, R. Schober and G. K. Karagiannidis, "Accelerating Distributed Optimization via Over-the-Air Computing," in IEEE Transactions on Communications, vol. 71, no. 9, pp. 5565-5579, Sept. 2023, doi: 10.1109/TCOMM.2023.3286915.
keywords: {Optimization;Convergence;Smart grids;Servers;Resource management;Linear programming;6G mobile communication;Over-the-air computing;non-orthogonal multiple access;primal-dual;distributed optimization;subgradient method;6G;large-scale optimization},

## License

- The code in this repository is licensed under the MIT License. See the LICENSE file for details.
- The paper and any associated figures or text are licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) License.
