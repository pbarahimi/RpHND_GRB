package RpHND_GRB.model;

public class Distance {

	public static double[][] get(double[][] input) {
		double[][] distance = new double[input.length][input.length];
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input.length; j++) {
				distance[i][j] = Math.sqrt((Math.pow(input[i][0] - input[j][0],
						2) + Math.pow(input[i][1] - input[j][1], 2)));
			}
		}
		return distance;
	}

	public static double[][] get(int[][] input) {
		double[][] distance = new double[input.length][input.length];
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input.length; j++) {
				distance[i][j] = Math.sqrt((Math.pow(input[i][0] - input[j][0],
						2) + Math.pow(input[i][1] - input[j][1], 2)));
			}
		}
		return distance;
	}

	public static double[][] get(float[][] input) {
		double[][] distance = new double[input.length][input.length];
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input.length; j++) {
				distance[i][j] = Math.sqrt((Math.pow(input[i][0] - input[j][0],
						2) + Math.pow(input[i][1] - input[j][1], 2)));
			}
		}
		return distance;
	}
}