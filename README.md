# 《数值代数》课程代码

教材：徐树方. 高立. 张平文 《数值线性代数（第二版）》

## 代码目录

- 第一章
  - 上三角形方程求解 `solveUpperTriangularEquation.m`
  - 下三角形方程求解 `solveLowerTriangularEquation.m`
  - 单位下三角形方程求解 `solveUnitLowerTriangularEquation.m`
  - LU分解 `getLU.m`
  - 列主元法LU分解 `getPLU.m`
  - Cholesky分解（LL分解） `cholesky.m`
  - 改进Cholesky分解（LDL分解） `choleskyImproved.m`
  - Gill-Murray修正Cholesky分解 `cholesky_GM.m`
  - 使用LU分解法解方程 `solveEquationWithLU.m`
  - 使用列主元LU分解法解方程 `solveEquationWithPLU.m`
  - 使用Cholesky分解法解方程 `solveEquationWithCholesky.m`
  - 使用改进Cholesky分解法解方程 `solveEquationWithCholeskyImproved.m`
  - 上机习题 `ex_1_1.m`、`ex_1_2_1.m`、`ex_1_2_2.m`
- 第二章
  - 优化法估计矩阵的一范数 `norm1Optimize.m`
  - 求矩阵的一范数条件数 `kappa.m`
  - 上机习题 `ex_2_1_1.m`
- 第三章
  - Householder变换 `householder.m`
  - Givens变换 `givens.m`
  - QR分解 `getQR.m`
  - 使用QR分解求解最小二乘问题 `solveLSwithQR.m`
  - 上机习题 `ex_3_1_1.m`、`ex_3_1_2.m`、`ex_3_1_3.m`、`ex_3_2.m`、`ex_3_3.m`
- 第四章
  - Jacobi迭代法 `jacobi.m`
  - Gauss-Seidel迭代法 `gauss_seidel.m`
  - 超松弛迭代法 `sor.m`
  - 上机习题 `ex_4_1.m`
- 第五章
  - 共轭梯度法 `CG.m`
  - 上机习题 `ex_5_2.m`、`ex_5_3.m`
- 第六章
  - 幂法求模最大特征值 `maxeig.m`
  - 用幂法求多项式模最大根 `maxroot.m`
  - 反幂法求特征向量 `inversePM.m`
  - 上Hessenberg分解 `hessenberg.m`
  - 双重步位移隐式QR迭代 `doubleQR.m`
  - 非对称实Schur分解 `realSchur.m`
  - 求非对称实矩阵的所有特征根 `eigen.m`
  - 求多项式的所有根 `allroot.m`
  - 上机习题 `ex_6_1.m`、`ex_6_2_2.m`、`ex_6_2_3.m`
- 第七章
  - 实对称矩阵的三对角化 `tridiag.m`
  - 带Wilkinson位移的隐式对称QR迭代 `wilkinsonQR.m`
  - 隐式对称QR迭代法求实对称矩阵的所有特征根和特征向量 `symmetricEigen.m`
  - 上机习题 `ex_7_1.m`
