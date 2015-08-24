import graph.Axis;
import graph.DataSet;
import graph.G2Dint;

import java.applet.Applet;
import java.awt.Button;
import java.awt.Color;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextArea;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Jama.Matrix;

/**
 * Jacobi.java, a class designed to run the Jacobi algorithm on a random matrix
 * and show the results in a graph.
 * 
 * @author Gene Hynson
 */
@SuppressWarnings("serial")
public class Jacobi extends Applet implements Runnable {
	/** The set for the actual results **/
	private DataSet data = new DataSet();
	/** The data set for the theoretical bound results **/
	private DataSet theorydata = new DataSet();
	/** The data set for the actual results when not sorting **/
	private DataSet noSort = new DataSet();
	/** The Off(A) numbers for each iteration **/
	private double[] offs;
	/** The natural log of 9/10 **/
	private final static double lognum = -0.10536051565782630122750098083931;
	private TextArea origDisp;
	private TextArea diagDisp;
	private TextArea noSortDisp;
	private G2Dint graph = new G2Dint();
	private Axis xaxis;
	private Axis yaxis;
	private Label stepsSort;
	private Label stepsNoSort;
	private final static String numStepsSort = "Steps to diagnolization when sorting: ";
	private final static String numStepsNoSort = "Steps to diagonlization when not sorting: ";
	private double stopAt = 10e-9;
	private int size = 5;
	private Matrix origMatrix;
	private Matrix diagMatrix;
	private boolean sorting = true;
	private Jacobi thisObject;
	
	public static void main(String[]args) {
		Jacobi j = new Jacobi();
		j.init();
		j.start();

	}
	/**
	 * Initialize the applet.
	 */
	@Override
	public void init() {
		super.init();
		setLayout(null);
		setSize(680, 460);
		setBackground(Color.WHITE);
		stepsSort = new Label();
		stepsSort.setBounds(320, 280, 300, 25);
		add(stepsSort);
		stepsNoSort = new Label();
		stepsNoSort.setBounds(320, 300, 300, 25);
		add(stepsNoSort);
		thisObject = this;
		makeMatrixBoxes();
		makeMatrix();
		
        new Thread(thisObject).start();            
		
		Panel button = new Panel();
        button.setBounds(100, 300, 500, 70);
        button.setLayout(null);
        add(button);
        
		Button repeat = new Button("Run Again With New Matrix");
        repeat.addActionListener(new ActionListener(){
           @Override
		public void actionPerformed(ActionEvent e) {
                    makeMatrix();
                    run();
           }
        });
        repeat.setBounds(275, 40, 160, 25);
        button.add(repeat);
		graph();

	}

	/**
	 * Makes the text boxes for the three matrices.
	 */
	private void makeMatrixBoxes() {
		Label orig = new Label("Original Matrix");
		orig.setBounds(5, 2, 100, 20);
		add(orig);
		origDisp = new TextArea();
		origDisp.setBounds(5, 25, 265, 115);
		origDisp.setEditable(false);
		add(origDisp);

		Label diag = new Label("Diagnolized Matrix");
		diag.setBounds(5, 155, 110, 20);
		add(diag);
		diagDisp = new TextArea();
		diagDisp.setEditable(false);
		diagDisp.setBounds(5, 178, 265, 115);
		add(diagDisp);

		Label nosort = new Label("Non-Sorted Matrix");
		nosort.setBounds(5, 308, 110, 20);
		add(nosort);
		noSortDisp = new TextArea();
		noSortDisp.setEditable(false);
		noSortDisp.setBounds(5, 331, 265, 115);
		add(noSortDisp);
	}

	/**
	 * Create the graph.
	 */
	public void graph() {

		xaxis = graph.createXAxis();
		xaxis.setTitleText("X");

		yaxis = graph.createYAxis();
		yaxis.setTitleText("xln(9/10)+ln(Off(A)) & ln(Off(B)");

		data.linestyle = 1;
		data.linecolor = Color.GREEN;
		data.legend(100, 194, "Results");
		data.legendColor(Color.black);

		theorydata.linestyle = 1;
		theorydata.linecolor = Color.RED;
		theorydata.legend(100, 210, "Theoretical Bound");
		theorydata.legendColor(Color.black);

		noSort.linestyle = 1;
		noSort.linecolor = Color.BLACK;
		noSort.legend(100, 226, "Non-Sorting");
		noSort.legendColor(Color.black);

		graph.attachDataSet(data);
		graph.attachDataSet(theorydata);

		xaxis.attachDataSet(data);
		xaxis.attachDataSet(theorydata);

		yaxis.attachDataSet(data);
		yaxis.attachDataSet(theorydata);

		graph.setDataBackground(Color.WHITE);
		graph.setBounds(260, 5, 425, 285);
		add(graph);
	}

	/**
	 * Draws the three graphs.
	 */
	public void plotStuff(DataSet data) {
		if (offs != null) {

			double[] array = new double[offs.length * 2];
			int i = 0;
			for (int j = 0; j < offs.length - 1; j++) {
				i += 2;
				array[i] = j;
				array[i + 1] = Math.log(offs[j]);
			}

			double[] points = new double[offs.length * 2];
			i = 0;
			for (int j = 0; j < offs.length - 1; j++) {
				i += 2;
				points[i] = j;
				points[i + 1] = (j * lognum) + Math.log(offs[0]);
			}

			data.deleteData();
			theorydata.deleteData();

			try {
				data.append(array, array.length / 2);
				theorydata.append(points, points.length / 2);
			} catch (Exception e) {
				e.printStackTrace();
			}

			xaxis.attachDataSet(data);
			xaxis.attachDataSet(theorydata);

			yaxis.attachDataSet(data);
			yaxis.attachDataSet(theorydata);

			graph.attachDataSet(data);
			graph.attachDataSet(theorydata);

			return;
		}
		graph.repaint();
	}

	/**
	 * Makes a random Matrix
	 */
	public void makeMatrix() {
		origMatrix = diagMatrix = Matrix.random(size, size);
		//symmetric matrix 
		for (int m = 0; m < size; m++) {
			for (int n = m + 1; n < size; n++) {
				origMatrix.set(n, m, origMatrix.get(m, n));
			}
		}
		graph.detachDataSets();
		graph.repaint();

		stepsSort.setText("");
		stepsNoSort.setText("");

		origDisp.setText(toString(origMatrix));
	}

	/**
	 * Returns the matrix as a string so that it can be displayed inside a text
	 * box.
	 */
	private String toString(Matrix matrix) {
		ByteArrayOutputStream byout = new ByteArrayOutputStream();
		PrintWriter out = new PrintWriter(byout);
		matrix.print(out, 10, 3);
		out.flush();

		String patString = "\n       ";
		String replaceString = "\n ";
		Pattern pat = Pattern.compile(patString);
		Matcher mat = pat.matcher(byout.toString());
		String matrixString = mat.replaceAll(replaceString);

		patString = "\n      -";
		replaceString = "\n-";
		pat = Pattern.compile(patString);
		mat = pat.matcher(matrixString);
		matrixString = mat.replaceAll(replaceString);

		matrixString = matrixString.replaceFirst("\n", "");

		return matrixString;
	}
	
	/**
	 * Runs the configuration.
	 */
	@Override
	public void run() {
		// run with sorting
		sorting = true;
		jacobi();
		diagDisp.setText(toString(diagMatrix));
		stepsSort.setText(numStepsSort + (offs.length - 1));
		plotStuff(data);

		// re-run without sorting
		sorting = false;
		jacobi();
		noSortDisp.setText(toString(diagMatrix));
		stepsNoSort.setText(numStepsNoSort + (offs.length - 1));
		plotStuff(noSort);

	}

	/**
	 * Performs the Jacobi on the matrix
	 */
	public void jacobi() {
		diagMatrix = origMatrix;
		List<Double> offlist = new ArrayList<Double>();
		int[] largestlocation = new int[] { 1, 0 };
		boolean keepGoing = true;
		while (keepGoing == true) {
			// Check Off(diagMatrix)
			// sum the squares of the off-diagnols
			double off = 0, largestoff = 0, temp;
			// calculate Off(diagMatrix) and find the largest off element
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					if (i != j) {
						temp = Math.abs(diagMatrix.get(i, j));
						if (sorting) {
							// check if the current entry is
							// larger than previously larger entry
							if (temp > largestoff) {
								largestoff = temp;
								// save location of entry
								largestlocation[0] = i;
								largestlocation[1] = j;
							}
						}
						off += Math.pow(temp, 2);
					}
				}
			}

			offlist.add(new Double(off));
			// check that off is <= 10^-9
			if (off <= stopAt) {
				keepGoing = false;
			}
			int x, y;
			// get the inner 2x2 matrix with the largest entry
			// or the next systematic matrix
			if (!sorting) {
				// last column, increment row, reset column
				if (largestlocation[0] == (size - 1)) {
					largestlocation[1] += 1;
					largestlocation[0] = largestlocation[1] + 1;
				} else {
					largestlocation[0] += 1;
				}
				// Once at the end of column, repeat
				if (largestlocation[1] == (size - 1)) {
					largestlocation = new int[] { 1, 0 };
				}
			}
			x = largestlocation[0];
			y = largestlocation[1];			
			Matrix k = new Matrix(2,2);
			k.set(1, 0, diagMatrix.get(x,y));
			k.set(0, 1, diagMatrix.get(x,y));
			k.set(1, 1, diagMatrix.get(y,y));
			k.set(0, 0, diagMatrix.get(x,x));			
			//k = findEigen(k); 
			k = k.eig().getV();

			// create the rotation matrix
			Matrix r = Matrix.identity(size, size);
			r.set(x, x, k.get(0, 0));
			r.set(x, y, k.get(0, 1));
			r.set(y, x, k.get(1, 0));
			r.set(y, y, k.get(1, 1));

			// produce the more diagnolized matrix
			diagMatrix = (r.transpose()).times(diagMatrix).times(r);
		}
		int length = offlist.size();
		offs = new double[length];
		for (int i = 0; i < length; i++) {
			offs[i] = offlist.get(i).doubleValue();
		}
	}
	//could not make this work
	private Matrix findEigen(Matrix q) {
		//finds eigenvalues
		double eig1 = ((q.get(1,1)+q.get(0, 0)) + Math.sqrt(Math.pow((q.get(1,1)-q.get(0,0)),2) + 4*Math.pow(q.get(0,1),2)))/2;
		//finds eigenvector
		Matrix i = Matrix.identity(2, 2);	
		Matrix k1 = q.minus(i.times(eig1));
		double[][] a1 = k1.getArray();
		double[][] b = new double[2][2];
		b[0][0] = a1[0][0];
		b[1][0] = a1[0][1];
		b[0][1] = -(a1[0][1]);
		b[1][1] = a1[0][0];		
		//normalizes eigenvectors
		double mag = Math.sqrt(Math.pow(a1[0][0],2) + Math.pow(a1[1][0],2));
		for (int j = 0; j < 2; j++) {
			for (int r = 0; r < 2; r++) {
				b[j][r] = b[j][r]/mag;
				System.out.println(b[j][r]);
				}	
			}
		q.equals(b);
		return q;
	}
}