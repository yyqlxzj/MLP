1 Introduction:
The purpose of this procedure is to calculate the maximum likelihood proton path based on the information of proton outgoing and incoming.
If necessary, please refer to the paper, named as A maximum likelihood proton path formalism for application in proton computed tomography, 
written by R. W. Schulte, for the specific form and derivation of the MLP formula.
Assuming that you have some problems during use, you can contact me by email,2107217959@qq.com.

2 Function:
func1,func2 and func3 are integrand functions in the covariance, whose parameters are fitted by GEANT4.
Input of function MLP includes the depth to be calculated and the position and angle at which the proton enters and exits.
Output of function MLP is a 2*2 matrix whose first row includes the positions of x and y respectively and second row includes the angles of
x and y respectively.You can define a new matrix in the main function to store the matrix calculated by the MLP

3 Suggestion:
You should pay attention to the defined domain of the func1,2,3 and u0,u1,u2 in log. Otherwise, you may get infinity, which will consume 
a lot of time to get a wrong result.
