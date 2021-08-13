# Secured-Constraints-Unit-Commitment-SCUC-model-of-Power-system  
Secured Constraints Unit Commitment model of Power system, including the model based on AC flow equation and DC flow equation(the AC model is relaxed by 'Second Order Cone Relaxation' to get a convex problem). 
Only the constraints in pre-contingency state are considered. 
There is a piecewise function expression for power generation costs.   
However, I'm sorry that the version is relative old, and it may be difficult to read or modify.   
电力系统安全约束机组组合问题模型，包括基于交流潮流方程和直流潮流方程的模型（交流潮流模型进行了二阶锥松弛以确保问题是凸的）。  
仅考虑意外前状态下的约束。  
有一个用于发电成本的分段函数表达式。  
但是，很抱歉该版本相对较旧，可能难以阅读或修改。  

Besides, note that the model is based on the Matlab, Yalmip, and the solver is Gurobi. It can be changed to other solvers, such as Cplex, by modifying the parameter 'gurobi' in sentence 'ops = settings('solver','gurobi''.   
此外，请注意，该模型基于Matlab，Yalmip，求解器为Gurobi。 通过修改句子'ops = settings('solver'，'gurobi'中的参数'gurobi'，可以将其更改为其他求解器，例如Cplex。  

If you have any idea on improving this model, please contact me.   
如果您有任何改进此模型的想法，请联系我.

注：matlab上传github，中文注释乱码的问题，暂时没找到解决办法，如果您有，可以告诉我。目前我的解决方法是，上传了一个.zip的压缩包，下载后解压应该不会乱码。
Note that: when I upload the .m files to github, its Chinese notes will be error codes, I didn't find any way to solve it, if you know, please tell me. Now my solution is that a .zip file is uploaded, upzip it seems to avoid the error codes problem.
