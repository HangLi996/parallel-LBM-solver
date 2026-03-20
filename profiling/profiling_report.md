# SIMD & Performance Profiling Report

## 1. Compiler Vectorization Analysis
**Status**: Mixed Results (Heavy fragmentation)

`vec_report.txt` confirms why:
-   **Memory Allocation in Loops**: `operator new` calls (from `std::vector` resizing/creation) are blocking vectorization.
-   **Control Flow**: complex branching inside `stream()`.

---

## 3. Likwid Profiling Results
**Group**: FLOPS_DP

| Metric | Value | Note |
| :--- | :--- | :--- |
| **DP [MFLOP/s]** | **36.24** | Extremely low. |
| **AVX DP [MFLOP/s]** | **4.35** | Only ~12% efficient. |
| **AVX512 DP [MFLOP/s]** | **0.00** | **NO AVX-512 usage.** |
| **Vectorization Ratio** | **3.87 %** | 🚨 **Critical Issue**. |

---

## 4. Recommendations
The code is currently **scalar-bound** due to memory allocations in the hot loop.

1.  **Refactor**: Move `std::vector` allocations out of `stream()` and `collide()`.
2.  **Simplify**: Use a branchless streaming (pull) method.
3.  **Manual Vectorization**: Only if refactoring fails.
