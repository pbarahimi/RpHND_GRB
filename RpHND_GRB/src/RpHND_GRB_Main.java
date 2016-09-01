import java.io.FileNotFoundException;

import RpHND_GRB.model.BinaryTree;
import RpHND_GRB.model.Distance;
import gurobi.*;

public class RpHND_GRB_Main {
	private static double alpha = 0.2;
	private static double[][] tmpFlows = MyArray.read("w.txt");
	private static double[][] coordinates = MyArray.read("coordinates.txt");
	private static double[][] fixedCosts = MyArray.read("fixedcharge.txt");
	private static double[][] distances = Distance.get(coordinates);//MyArray.read("CAB10_Distance.txt");
	private static int nVar = tmpFlows.length;
	private static double[][] flows = new double[nVar][nVar];
	private static int D = 1; // Maximum number of simultaneous disruptions
	private static int R = (int) Math.pow(2, D + 1) - 2; // Largest index in the
															// full binary tree
	private static double q = 0.2;
	private static int P = 3; // number of hubs to be located
	private static int M = nVar * R; // the big M

	public static void main(String[] args) throws FileNotFoundException {
		
		//Filling in the flows matrix assymetrically 
		for (int i = 0; i < nVar; i++) {
			for (int j = 0; j < nVar; j++) {
				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
			}
		}
		
		try {
			GRBEnv env = new GRBEnv(null);
			GRBModel model = new GRBModel(env);

			// Create variables
			GRBVar[][][][][] x = new GRBVar[nVar][nVar][nVar][nVar][R + 1];
			GRBVar[] y = new GRBVar[nVar];

			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int m = 0; m < nVar; m++) {					
							String varName = "x" + i + "_" + k + "_" + m + "_"
									+ j + "_0";
							x[i][k][m][j][0] = model.addVar(0.0, 1.0, 0.0,
									GRB.BINARY, varName);
						}
					}
				}
			}		

			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int m = 0; m < nVar; m++) {
							for (int r = 1; r <= R; r++){
								String varName = "x" + i + "_" + k + "_" + m
										+ "_" + j + "_" + r;
								x[i][k][m][j][r] = model.addVar(0.0, 1.0, 0.0,
										GRB.BINARY, varName);
							}
						}
					}
				}
			}

			for (int i = 0; i < nVar; i++) {
				y[i] = model.addVar(0, 1.0, 0, GRB.BINARY, "y" + i);
			}

			// Integrate new variables
			model.update();

			// Set objective function
			GRBLinExpr expr = new GRBLinExpr();
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int k = 0; k < nVar; k++) {
						for (int m = 0; m < nVar; m++) {
							double CoEf = -1 * flows[i][j] * Cikmj(i, k, m, j) * (1 - Q(i,k,m,j));
							expr.addTerm(CoEf, x[i][k][m][j][0]);
							if (i==3 && k==2&&m==2&&j==4)
								System.out.println(Cikmj(i, k, m, j) + "-" + (1 - Q(i,k,m,j)));
						}
					}
				}
			}

			for (int r = 1; r <= R; r++) {
				for (int i = 0; i < nVar; i++) {
					for (int j = i + 1; j < nVar; j++) {
						for (int k = 0; k < nVar; k++) {
							for (int m = 0; m < nVar; m++) {
								double CoEf = -1 * flows[i][j] * Cikmj(i, k, m, j) * Math.pow(q, Math.floor(Math.log(r+1)/Math.log(2)));
								expr.addTerm(CoEf, x[i][k][m][j][r]);
							}
						}
					}
				}
			}
			
			for (int i = 0; i < nVar; i++) {
				expr.addTerm(-1* fixedCosts[i][0], y[i]);
			}

			model.setObjective(expr, GRB.MAXIMIZE);

			// Adding constraints
			// Constraint 2
			GRBLinExpr con2 = new GRBLinExpr();
			for (int i = 0; i < nVar; i++) {
				con2.addTerm(1, y[i]);
			}
			model.addConstr(con2, GRB.EQUAL, P, "u2");

			// Constraint 3
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {

					GRBLinExpr con3 = new GRBLinExpr();
					for (int k = 0; k < nVar; k++) {
						for (int m = 0; m < nVar; m++) {
							con3.addTerm(1, x[i][k][m][j][0]);
						}
					}
					model.addConstr(con3, GRB.EQUAL, 1, "u3_" + i + "_" + j);
				}
			}

			// Constraint 4
			for (int r = 0; r <= R; r++) {
				for (int i = 0; i < nVar; i++) {
					for (int j = i + 1; j < nVar; j++) {
						for (int k = 0; k < nVar; k++) {

							GRBLinExpr con4 = new GRBLinExpr();
							for (int m = 0; m < nVar; m++) {
								con4.addTerm(1, x[i][k][m][j][r]);
							}
							con4.addTerm(-1, y[k]);
							model.addConstr(con4, GRB.LESS_EQUAL, 0, "u4_" + i + "_" + j + "_" + k + "_" + r);
						}
					}
				}
			}

			// Constraint 5
			for (int r = 0; r <= R; r++) {
				for (int i = 0; i < nVar; i++) {
					for (int j = i + 1; j < nVar; j++) {
						for (int m = 0; m < nVar; m++) {

							GRBLinExpr con5 = new GRBLinExpr();
							for (int k = 0; k < nVar; k++) {
								con5.addTerm(1, x[i][k][m][j][r]);
							}
							con5.addTerm(-1, y[m]);
							model.addConstr(con5, GRB.LESS_EQUAL, 0, "u5_" + i + "_" + j + "_" + m + "_" + r);
						}
					}
				}
			}

			// Constraint 6
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {

						GRBLinExpr con6 = new GRBLinExpr();
						for (int k = 0; k < nVar; k++) {
							for (int m = 0; m < nVar; m++) {
								if (k != i && m != i) {
									con6.addTerm(1, x[i][k][m][j][r]);
								}
							}
						}
						con6.addTerm(M, y[i]);
						model.addConstr(con6, GRB.LESS_EQUAL, M, "u6_" + i + "_" + j + "_" + r);
					}
				}
			}

			// Constraint 7
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= R; r++) {

						GRBLinExpr con7 = new GRBLinExpr();
						for (int k = 0; k < nVar; k++) {
							for (int m = 0; m < nVar; m++) {
								if (k != j && m != j) {
									con7.addTerm(1, x[i][k][m][j][r]);
								}
							}
						}
						con7.addTerm(M, y[j]);
						model.addConstr(con7, GRB.LESS_EQUAL, M, "u7_" + i + "_" + j + "_" + r);
					}
				}
			}

			// Constraint 8
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= Math.pow(2, D) - 2; r++) { // The
																	// leaf-nodes
																	// are not
																	// to be
																	// considered
																	// in this
																	// constraint.

						GRBLinExpr con8 = new GRBLinExpr();
						for (int k = 0; k < nVar; k++) {
							for (int m = 0; m < nVar; m++) {
								if (k != i && k != j) {
									con8.addTerm(1, x[i][k][m][j][r]);
								}
							}
						}

						for (int k = 0; k < nVar; k++) {
							for (int m = 0; m < nVar; m++) {
								con8.addTerm(-1, x[i][k][m][j][2 * r + 1]); // left
																			// child
																			// node

							}
						}
						model.addConstr(con8, GRB.LESS_EQUAL, 0, "u8_" + i + "_" + j + "_" + r);
					}
				}
			}

			// Constraint 9
			for (int i = 0; i < nVar; i++) {
				for (int j = i + 1; j < nVar; j++) {
					for (int r = 0; r <= Math.pow(2, D) - 2; r++) { // The
																	// leaf-nodes
																	// are not
																	// to be
																	// considered
																	// in this
																	// constraint.

						GRBLinExpr con9 = new GRBLinExpr();
						for (int k = 0; k < nVar; k++) {
							for (int m = 0; m < nVar; m++) {
								if (m!=i && m!=j && m!=k){
									con9.addTerm(1, x[i][k][m][j][r]);
								}
							}
						}
						
						for (int k = 0; k < nVar; k++) {
							for (int m = 0; m < nVar; m++) {
								con9.addTerm(-1, x[i][k][m][j][2 * r + 2]); // right
																			// child
																			// node
							}
						}
						model.addConstr(con9, GRB.LESS_EQUAL, 0, "u9_" + i + "_" + j + "_" + r);
					}
				}
			}

			// Constraint 10
			for (int i=0;i<nVar;i++){
				for (int j=i+1;j<nVar;j++){
					for (int k=0;k<nVar;k++){
						for (int r=0;r<=Math.pow(2, D) - 2;r++){
							GRBLinExpr con10 = new GRBLinExpr();
							for (int s:BinaryTree.getLeftChildren(r, D)){
								for (int m=0; m<nVar; m++){
									con10.addTerm(1, x[i][k][m][j][s]);
									con10.addTerm(1, x[i][m][k][j][s]);
								}
							}
							for (int m=0; m<nVar; m++){
								con10.addTerm(M , x[i][k][m][j][r]);
							}
							model.addConstr(con10, GRB.LESS_EQUAL, M, "u10_" + i + "_" + j + "_" + k + "_" + r);
						}
					}
				}
			}

			// Constraint 11
			for (int i=0;i<nVar;i++){
				for (int j=i+1;j<nVar;j++){
					for (int m=0;m<nVar;m++){
						for (int r=0;r<=Math.pow(2, D) - 2;r++){

							GRBLinExpr con11 = new GRBLinExpr();
							for (int s : BinaryTree.getRightChildren(r, D)) {
								for (int k = 0; k < nVar; k++) {
									con11.addTerm(1, x[i][k][m][j][s]);
									con11.addTerm(1, x[i][m][k][j][s]);
								}
							}
							for (int k = 0; k < nVar; k++) {
								con11.addTerm(M, x[i][k][m][j][r]);
							}
							model.addConstr(con11, GRB.LESS_EQUAL, M, "u11_" + i + "_" + j + "_" + m + "_" + r);
						}
					}
				}
			}
			
			/*GRBLinExpr expr0 = new GRBLinExpr();
			expr0.addTerm(1, y[1]);
			model.addConstr(expr0, GRB.EQUAL, 1, null);
			
			expr0 = new GRBLinExpr();
			expr0.addTerm(1, y[2]);
			model.addConstr(expr0, GRB.EQUAL, 1, null);
			
			expr0 = new GRBLinExpr();
			expr0.addTerm(1, y[6]);
			model.addConstr(expr0, GRB.EQUAL, 1, null);*/
			
			
			// Optimize model
			double startTime = System.currentTimeMillis();
			model.optimize();
			System.out.println("Elapsed time: "  + (System.currentTimeMillis() - startTime));
						
			// Printing solution to a file
//			File file = new File("D:/primalResults.txt");
//			PrintWriter out = new PrintWriter(file);
			
				
			
			/*GRBVar[] vars = model.getVars();
			for (GRBVar var: vars){
				out.println(var.get(GRB.DoubleAttr.Obj) + " " + var.get(GRB.StringAttr.VarName) + " " + var.get(GRB.DoubleAttr.X));
			}*/
			
			/*GRBConstr[] constrs = model.getConstrs();
			for (GRBConstr constr: constrs){
				System.out.println(constr.get(GRB.StringAttr.ConstrName) + " " + constr.get(GRB.DoubleAttr.Pi));
			}*/
			
//			out.close();
			for (GRBVar var : model.getVars())
				System.out.println(var.get(GRB.StringAttr.VarName) + ":" + var.get(GRB.DoubleAttr.Obj));
			
			System.out.println("Number of variables: " + model.get(GRB.IntAttr.NumVars));
			System.out.println("Number of constraints: " + model.get(GRB.IntAttr.NumConstrs));
			printSol(model);
			
			//Results 
			System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));
			GRBVar var = model.getVarByName("x0_3_3_4_2");
			System.out.println(var.get(GRB.StringAttr.VarName) + ": " + var.get(GRB.DoubleAttr.Obj));
			// Dispose of model and environment
		    model.dispose();
		    env.dispose();

		    
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". "
					+ e.getMessage());
		}

	}

	
	/**
	 * prints the solution
	 * @param model
	 * @throws GRBException
	 */
	public static void printSol(GRBModel model) throws GRBException{
		for (GRBVar var : model.getVars()) 
			if (var.get(GRB.DoubleAttr.X)>0 /*&& var.get(GRB.StringAttr.VarName).contains("y")*/) 
				System.out.println(var.get(GRB.StringAttr.VarName) + " : " + var.get(GRB.DoubleAttr.X) + " : " + var.get(GRB.DoubleAttr.Obj));
	}
	
	/**
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @param m
	 * @return operating probability of a route
	 */
	public static double Q(int i, int k, int m, int j) {
		double result = 0;
		if (k!=i && j!=m)
				result = q+q(k,m);	
		else if (i==k && m!=j)
			result = q(i,m);
		else if (i!=k && m==j)
			result = q(j,m);
		else if (i==k && j==m)
			result = 0;
		else 
			System.out.println("Not include in the Q(i,k,m,j)!");
		return result;
	}

	/**
	 * Cikmj
	 * 
	 */
	private static double Cikmj(int i, int k, int m, int j) {
		double cost = distances[i][k] + (1 - alpha) * distances[k][m]
				+ distances[m][j];
		/*
		 * double cost = collCost * distances[i][k] + transCost *
		 * distances[k][m] + distCost * distances[m][j];
		 */
		return cost;
	}

	/**
	 * q(k,m)
	 */
	private static double q(int k, int m) {
		if (k == m)
			return 0;
		else
			return q;
	}
}