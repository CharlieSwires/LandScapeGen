/*
Author CEL Swires all rights reserved.

Non commercial use permitted.

Copyright 2010 CEL Swires
*/
import java.io.*;
import java.util.ArrayList;
import java.util.Random;
import java.lang.Math;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.awt.*;

import javax.imageio.ImageIO;
import javax.swing.*;

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.awt.event.*;
import java.awt.image.BufferedImage;

public class LandscapeGen extends JFrame{

	double sceneArray[][];

	static final double roughness = 0.3;
	Random r;
	double fa[];
	final double scale = (double)1000.0;


		JButton render, quit;
		JLabel JLseed,JLseed2, JLheightScale, JLh;
		JTextField JTFseed,JTFSeed2, JTFheightScale, JTFh;
		WindAdapter windListener;
		RenderCanvas ro;
		Container contentPane;
		long seed;long seed2; 
		String outFile;String outFile2;
		String inOutFile; 
		String heightFile; 
		int edgeLength; String landType;

		public LandscapeGen(String title, long seed, String outFile, int edgeLength, String landType, long seed2, String normalImage, String objFile, String heightFile){
			setTitle(title);
			this.seed = seed;this.seed2 = seed2; this.inOutFile = objFile;this.heightFile = heightFile;this.landType =landType;
			this.outFile2 = normalImage;this.outFile = outFile; this.edgeLength = edgeLength;
			contentPane = getContentPane();
			contentPane.setLayout(new BorderLayout());
			JPanel buttonPanel = new JPanel();
			buttonPanel.setLayout(new FlowLayout());
			ButtonListener listener = new ButtonListener();
			windListener = new WindAdapter();
			this.addWindowListener((WindowListener)windListener);
			render = new JButton("Render");
			render.addActionListener(listener);
			buttonPanel.add(render);
			quit = new JButton("Quit");
			quit.addActionListener(listener);
			buttonPanel.add(quit);
			contentPane.add(buttonPanel, BorderLayout.SOUTH);
			JPanel dimensionsPanel = new JPanel();
			dimensionsPanel.setLayout(new FlowLayout());
			JPanel seedPanel = new JPanel();
			seedPanel.setLayout(new FlowLayout());
			JLseed = new JLabel("seed:");
			JTFseed = new JTextField(""+seed);
			seedPanel.add(JLseed);
			seedPanel.add(JTFseed);
			JLseed2 = new JLabel("seed craters:");
			JTFSeed2 = new JTextField(""+seed2);
			seedPanel.add(JLseed2);
			seedPanel.add(JTFSeed2);
			dimensionsPanel.add(seedPanel);
			JPanel heightScalePanel = new JPanel();
			heightScalePanel.setLayout(new FlowLayout());
			JLheightScale = new JLabel("height scale:");
			JTFheightScale = new JTextField("1.0");
			heightScalePanel.add(JLheightScale);
			heightScalePanel.add(JTFheightScale);
			dimensionsPanel.add(heightScalePanel);
			JPanel hPanel = new JPanel();
			hPanel.setLayout(new FlowLayout());
			JLh = new JLabel("h:");
			JTFh = new JTextField("1.0");
			hPanel.add(JLh);
			hPanel.add(JTFh);
			dimensionsPanel.add(hPanel);
			contentPane.add(dimensionsPanel, BorderLayout.NORTH);
			ro = new RenderCanvas(seed, outFile, edgeLength);
			contentPane.add(ro, BorderLayout.CENTER);
			setSize(800,800);
			setVisible(true);
		}
	class ButtonListener implements ActionListener{
		public void actionPerformed(ActionEvent e){
			if (e.getSource() == render){
				contentPane.remove(ro);
				ro = new RenderCanvas(seed, outFile, edgeLength);
				contentPane.add(ro, BorderLayout.CENTER);
				setVisible(true);
				repaint();

			}
			if (e.getSource() == quit){
				System.exit(0);
			}
		}
	}

        class WindAdapter extends WindowAdapter{

                public void windowClosing(WindowEvent e){
 			System.exit(0);
                }

        }

/*
 Modified for Java by CEL Swires.

 Written by: Paul E. Martz

 Copyright 1997 by Paul E. Martz, all right reserved

 Non-commercial use by individuals is permitted.

*/

/*
 * randNum - Return a random doubleing point number such that
 *      (min <= return-value <= max)
 * 32,767 values are possible for any given range.
 */
double randnum (double min, double max)
{
    double	x;

	x = r.nextFloat();
    return (x * (max - min) + min);
}


/*
 * fractRand is a useful interface to randnum.
 */
double fractRand(double v){
	return(randnum (-v, v));
}



/*
 * avgEndpoints - Given the i location and a stride to the data
 * values, return the average those data values. "i" can be thought of
 * as the data value in the center of two line endpoints. We use
 * "stride" to get the values of the endpoints. Averaging them yields
 * the midpoint of the line.
 *
 * Called by fill1DFractArray.
 */
static double avgEndpoints (int i, int stride, double fa[])
{
    return ((double) (fa[i-stride] +
		     fa[i+stride]) * .5f);
}

/*
 * avgDiamondVals - Given the i,j location as the center of a diamond,
 * average the data values at the four corners of the diamond and
 * return it. "Stride" represents the distance from the diamond center
 * to a diamond corner.
 *
 * Called by fill2DFractArray.
 */
static double avgDiamondVals (int i, int j, int stride,
			     int size, int subSize, double fa[])
{
    /* In this diagram, our input stride is 1, the i,j location is
       indicated by "X", and the four value we want to average are
       "*"s:
           .   *   .

           *   X   *

           .   *   .
       */

    /* In order to support tiled surfaces which meet seamless at the
       edges (that is, they "wrap"), We need to be careful how we
       calculate averages when the i,j diamond center lies on an edge
       of the array. The first four 'if' clauses handle these
       cases. The final 'else' clause handles the general case (in
       which i,j is not on an edge).
     */
    if (i == 0)
	return ((double) (fa[(i*size) + j-stride] +
			 fa[(i*size) + j+stride] +
			 fa[((subSize-stride)*size) + j] +
			 fa[((i+stride)*size) + j]) * .25f);
    else if (i == size-1)
	return ((double) (fa[(i*size) + j-stride] +
			 fa[(i*size) + j+stride] +
			 fa[((i-stride)*size) + j] +
			 fa[((0+stride)*size) + j]) * .25f);
    else if (j == 0)
	return ((double) (fa[((i-stride)*size) + j] +
			 fa[((i+stride)*size) + j] +
			 fa[(i*size) + j+stride] +
			 fa[(i*size) + subSize-stride]) * .25f);
    else if (j == size-1)
	return ((double) (fa[((i-stride)*size) + j] +
			 fa[((i+stride)*size) + j] +
			 fa[(i*size) + j-stride] +
			 fa[(i*size) + 0+stride]) * .25f);
    else
	return ((double) (fa[((i-stride)*size) + j] +
			 fa[((i+stride)*size) + j] +
			 fa[(i*size) + j-stride] +
			 fa[(i*size) + j+stride]) * .25f);
}


/*
 * avgSquareVals - Given the i,j location as the center of a square,
 * average the data values at the four corners of the square and return
 * it. "Stride" represents half the length of one side of the square.
 *
 * Called by fill2DFractArray.
 */
static double avgSquareVals (int i, int j, int stride, int size, double fa[])
{
    /* In this diagram, our input stride is 1, the i,j location is
       indicated by "*", and the four value we want to average are
       "X"s:
           X   .   X

           .   *   .

           X   .   X
       */
    return ((double) (fa[((i-stride)*size) + j-stride] +
		     fa[((i-stride)*size) + j+stride] +
		     fa[((i+stride)*size) + j-stride] +
		     fa[((i+stride)*size) + j+stride]) * .25f);
}


/*
 * powerOf2 - Returns 1 if size is a power of 2. Returns 0 if size is
 * not a power of 2, or is zero.
 */
boolean powerOf2 (int size)
{
    int i, bitcount = 0;

    /* Note this code assumes that (sizeof(int)*8) will yield the
       number of bits in an int. Should be portable to most
       platforms. */
    for (i=0; i<4*8; i++)
	if ((size & (1<<i))>0)
	    bitcount++;
    if (bitcount == 1)
	/* One bit. Must be a power of 2. */
	return (true);
    else
	/* either size==0, or size not a power of 2. Sorry, Charlie. */
	return (false);
}


/*
 * fill2DFractArray - Use the diamond-square algorithm to tessalate a
 * grid of double values into a fractal height map.
 */
void fill2DFractArray (double fa[], int size,
		       long seedValue, double heightScale, double h)
{
    int	i, j;
    int	stride;
    boolean	oddline;
    int subSize;
	double ratio, scale;

    if (!powerOf2(size) || (size==1)) {
	/* We can't tesselate the array if it is not a power of 2. */
	return;
    }

    /* subSize is the dimension of the array in terms of connected line
       segments, while size is the dimension in terms of number of
       vertices. */
    subSize = size;
    size++;

    /* initialize random number generator */
    r.setSeed (seedValue);


	/* Set up our roughness constants.
	   Random numbers are always generated in the range 0.0 to 1.0.
	   'scale' is multiplied by the randum number.
	   'ratio' is multiplied by 'scale' after each iteration
	   to effectively reduce the randum number range.
	   */
	ratio = (double) Math.pow (2.0,-h);
	scale = heightScale * ratio;

    /* Seed the first four values. For example, in a 4x4 array, we
       would initialize the data points indicated by '*':

           *   .   .   .   *

           .   .   .   .   .

           .   .   .   .   .

           .   .   .   .   .

           *   .   .   .   *

       In terms of the "diamond-square" algorithm, this gives us
       "squares".

       We want the four corners of the array to have the same
       point. This will allow us to tile the arrays next to each other
       such that they join seemlessly. */

    stride = subSize / 2;
    fa[(0*size)+0] =
      fa[(subSize*size)+0] =
        fa[(subSize*size)+subSize] =
          fa[(0*size)+subSize] = 0.f;


    /* Now we add ever-increasing detail based on the "diamond" seeded
       values. We loop over stride, which gets cut in half at the
       bottom of the loop. Since it's an int, eventually division by 2
       will produce a zero result, terminating the loop. */
    while ( stride > 0) {
		/* Take the existing "square" data and produce "diamond"
		   data. On the first pass through with a 4x4 matrix, the
		   existing data is shown as "X"s, and we need to generate the
	       "*" now:

               X   .   .   .   X

               .   .   .   .   .

               .   .   *   .   .

               .   .   .   .   .

               X   .   .   .   X

	      It doesn't look like diamonds. What it actually is, for the
	      first pass, is the corners of four diamonds meeting at the
	      center of the array. */
		for (i=stride; i<subSize; i+=stride) {
			for (j=stride; j<subSize; j+=stride) {
				fa[(i * size) + j] =
					scale * fractRand (.5f) +
					avgSquareVals (i, j, stride, size, fa);
				j += stride;
			}
			i += stride;
		}

		/* Take the existing "diamond" data and make it into
	       "squares". Back to our 4X4 example: The first time we
	       encounter this code, the existing values are represented by
	       "X"s, and the values we want to generate here are "*"s:

               X   .   *   .   X

               .   .   .   .   .

               *   .   X   .   *

               .   .   .   .   .

               X   .   *   .   X

	       i and j represent our (x,y) position in the array. The
	       first value we want to generate is at (i=2,j=0), and we use
	       "oddline" and "stride" to increment j to the desired value.
	       */
		oddline = false;
		for (i=0; i<subSize; i+=stride) {
			if(oddline == false){
				oddline = true;
			}else{
				oddline = false;
			}
			for (j=0; j<subSize; j+=stride) {
				if ((oddline) && j==0) j+=stride;

				/* i and j are setup. Call avgDiamondVals with the
				   current position. It will return the average of the
				   surrounding diamond data points. */
				fa[(i * size) + j] =
					scale * fractRand (.5f) +
					avgDiamondVals (i, j, stride, size, subSize, fa);

				/* To wrap edges seamlessly, copy edge values around
				   to other side of array */
				if (i==0)
					fa[(subSize*size) + j] =
						fa[(i * size) + j];
				if (j==0)
					fa[(i*size) + subSize] =
						fa[(i * size) + j];

				j+=stride;
			}
		}

		/* reduce random number range. */
		scale *= ratio;
		stride >>= 1;
    }

}

static double [] genNormal (double x1, double y1, double z1,
		       double x2, double y2, double z2,
		       double x3, double y3, double z3)
{
    double	len;
    double	v1x, v1y, v1z;
    double	v2x, v2y, v2z;
    double   normal[];

    normal = new double[3];


    v1x = x2 - x1;
    v1y = y2 - y1;
    v1z = z2 - z1;
    v2x = x3 - x1;
    v2y = y3 - y1;
    v2z = z3 - z1;

    normal[0] = v1y*v2z - v1z*v2y;
    normal[1] = v1z*v2x - v1x*v2z;
    normal[2] = v1x*v2y - v1y*v2x;

    len = (double) Math.sqrt ((double)(normal[0]*normal[0] + normal[1]*normal[1] +
			normal[2]*normal[2]));

    normal[0] /= len;
    normal[1] /= len;
    normal[2] /= len;
    return (normal);
}



void draw3DTriangle (double x1, double y1, double z1,
		     double x2, double y2, double z2,
		     double x3, double y3, double z3, 
		     PrintWriter out, double normal[], Dimension d, 
		     Graphics g, 
		     int j, int i, BufferedImage img, int size)
{

	double lightSource[] = new double[3];
	lightSource[0] = 0.0;
	lightSource[1] = 1.0;
	lightSource[2] = 0.0;
	double shade = normal[0]*lightSource[0]+normal[1]*lightSource[1]+normal[2]*lightSource[2];
	shade = (shade > 1.0) ? 1.0 : shade;
	shade = (shade < 0.0) ? 0.0 : shade;
	Color col;
//		if (y1 > 0.06){
//			col = Color.WHITE;
//			if (shade < 0.9 && y1 < 0.075)
//				col = Color.GRAY;
//		}
//		else if (y1 <= 0.06 && y1 > 0.0005){
//			col = Color.GRAY;
//	//		if (shade < 0.8)
//	//			col = Color.GRAY;
//		}
//		else if (y1 <= 0.0005 && y1 > -0.01){
//			col = Color.YELLOW;
//		}
//		else col = Color.BLUE;
	//col = Color.WHITE;
	col = Color.GRAY;
	Color c = new Color((float)(shade * col.getRed()/255.0),
			(float)(shade * col.getGreen()/255.0),
			(float)(shade* col.getBlue()/255.0));
	int xPoints[];
	int yPoints[];
	xPoints = new int[3];
	yPoints = new int[3];

	xPoints[0] = (int)((x1+1.0)*d.width/2.0);
	xPoints[1] = (int)((x2+1.0)*d.width/2.0);
	xPoints[2] = (int)((x3+1.0)*d.width/2.0);
	yPoints[0] = (int)((z1+1.0)*d.height/2.0);
	yPoints[1] = (int)((z2+1.0)*d.height/2.0);
	yPoints[2] = (int)((z3+1.0)*d.height/2.0);
	if (img != null){
	g.setColor(c);
	g.fillPolygon(xPoints, yPoints, 3);
	//shade = 1.0;
	c = new Color((float)(shade * col.getRed()/255.0),
			(float)(shade * col.getGreen()/255.0),
			(float)(shade* col.getBlue()/255.0));
	img.setRGB(j, i,c.getRGB());
	}//img.setRGB(size+(size-j)-1, i, c.getRGB());

	if(out!=null)out.println("tr "+ x1 * scale + ", " + z1 * scale + ", " + -y1 * scale + ", " +
	x2 * scale + ", " + z2 * scale + ", " + -y2* scale + ", " +
	x3 * scale + ", " + z3 * scale + ", " + -y3* scale + ";");

}
void bounding (double x1, double y1, double z1,
	     double x2, double y2, double z2,
	     double x3, double y3, double z3, 
	     double x4, double y4, double z4, 
	     double x5, double y5, double z5, 
	     PrintWriter out)
{
	if(out!=null)out.println("bounding "+ x1 * scale + ", " + z1 * scale + ", " + -y1 * scale + ", " +
	x2 * scale + ", " + z2 * scale + ", " + -y2* scale + ", " +
	x3 * scale + ", " + z3 * scale + ", " + -y3* scale + ", " +
	x4 * scale + ", " + z4 * scale + ", " + -y4* scale + ", " +
	x5 * scale + ", " + z5 * scale + ", " + -y5* scale +
	";");
}

/*
 * draw2DFractArrayAsTriangles - Draws the height map as a set of
 * triangular facets with a corresponding facet normal.
 *
 * This is a simplified routine intended for getting things up and
 * running quickly, and as a demonstration of how to walk through the
 * array.
 *
 * To use this routine, you MUST define your own "draw3DTriangle"
 * function according to the extern definition below. It takes 12
 * (ugh!) parameters: the first 9 are the X, Y, and Z world
 * coordinates of the three vertices, and the last three parameters
 * are the X, Y, and Z cartesian components of the unit-length facet
 * normal.
 *
 * X and Z coordinates will be distributed evenly over a grid from
 * -1.0 to 1.0 along both axes. Corresponding Y coordinates will be
 * extracted from the fract array. Normals will be gererated favoring
 * the "UP" or positive Y direction, and vertex order will use the
 * h hand rule to match the normal. ("h hand rule": When
 * viewed from the direction the normal is pointing, the vertex
 * ordering will be counter-clockwise.)
 */
void draw2DFractArrayAsTriangles (long seed2, double fa[], int size, PrintWriter out, Dimension d, Graphics g) throws IOException
{
    int	i, j;
    double	x, z, inc;
    final int	subSize = size;
	BufferedImage img = new BufferedImage(size, size,
			BufferedImage.TYPE_INT_RGB);
    size++;

    //set floor
    inc = 2.f / subSize;
	for (i = 0;i < subSize;i++){
		for (j = 0; j <subSize;j++){
			x = -1.0 + j*inc;
			z = -1.0 + i*inc;
			double temp = fa[(i*size)+j];
		    if (temp < 0.0){
				fa[(i*size)+j] = temp = 0.0;
			}    
		}
	}
	if (this.landType.equals("moon")){
	//generate slice
	BufferedImage craterCrossImg = ImageIO.read(new File("craterCrossSection.png"));
	double craterCrossArray[] = new double[subSize];
	for (i=0;i<subSize;i++){
		for (j = 0; j <craterCrossImg.getHeight();j++){
			if ((craterCrossImg.getRGB((int)(i*(double)craterCrossImg.getWidth()/subSize), j)& 0xff)>Color.GRAY.getBlue()){
				craterCrossArray[i]=craterCrossImg.getHeight()-(double)j;
				System.out.println(craterCrossImg.getHeight()-(double)j);
				break;
			}
		}
	}
	//normalise
	double zeroValue = craterCrossArray[subSize-1];
	for (i=0;i<subSize;i++){
		craterCrossArray[i] = craterCrossArray[i]-zeroValue;
	}
	// Collect data.
	final WeightedObservedPoints obs = new WeightedObservedPoints();
	for (i=0;i<subSize;i++){
		obs.add(i,craterCrossArray[i]);
	}

	// Instantiate a third-degree polynomial fitter.
	final PolynomialCurveFitter fitter = PolynomialCurveFitter.create(40);

	// Retrieve fitted parameters (coefficients of the polynomial function).
	final double[] coeff = fitter.fit(obs.toList());
	Random rand = new Random(seed2);
	//add craters
	for (int cratersi = 0; cratersi < 200; cratersi++){
		double centrex = rand.nextDouble() * subSize;
		double centrey = rand.nextDouble() * subSize;
		double r = Math.pow(10.0,((Math.log10(rand.nextDouble()*Math.pow(10.0, 4.0))-5.46)/-1.82));
		double centreHeight = fa[(int)((centrex*size)+centrey)];
		for (i = 0; i < subSize;i++){
			for (j=0;j<subSize;j++){
				double relx = i - centrex;
				double rely = j - centrey;
				if (rely * rely + relx * relx < r * r){
					fa[(i*size)+j] = fa[(i*size)+j]+r*0.01*craterCrossArray(coeff,(double)(subSize*((Math.sqrt((rely * rely + relx * relx)/(r*r))/2.0)+0.5)))/craterCrossImg.getWidth();
				}
			}
		}
	}
	//add curve
	for (i = 0;i < subSize;i++){
		for (j = 0; j <subSize;j++){
			x = -1.0 + j*inc;
			z = -1.0 + i*inc;
			double r = 1738.1 / 400.0;
			double t = fa[(i*size)+j];
			fa[(i*size)+j] = Math.sqrt((t+r)*(t+r)-
							x*x - z*z) - r;
		}
	}
	}
    divide (fa, subSize, 0,0,subSize,out,d,g,img);
    loadObjFileWithVertsAndNormals(this.inOutFile,fa,subSize);
    loadObjFileWithHeights(this.heightFile,fa,subSize);
//	System.out.println("bounding ;");
//	double xx = -5.0;
//	double zz = -5.0;
//	double xxw = 5.0;
//	double zzw = 5.0;
//	double min = 0.0001;
//	double max = -0.0001;
//	bounding ( xx,  min,  zz,
//		      xxw,  min,  zz,
//		      xxw,  min,  zzw, 
//		      xx,  min,  zzw, 
//		      xx,  max,  zz, 
//		     out);
//
//    double normal[] = genNormal (-5.0, 0.0, -5.0,
//		    -5.0, 0.0, 5.0,
//		    5.0, 0.0, -5.0);
//    draw3DTriangle (-5.0, 0.0, -5.0,
//		    -5.0, 0.0, 5.0,
//		    5.0, 0.0, -5.0, 
//		    out, normal, d, g,0,0,null,subSize);
//
//    double normal2[] = genNormal (-5.0, 0f, 5.0,
//		    5.0, 0.0, 5.0,
//		    5.0, 0.0, -5.0);
//
//    draw3DTriangle (-5.0, 0.0, 5.0,
//		    5.0, 0.0, 5.0,
//		    5.0, 0.0, -5.0, 
//		    out, normal2, d, g,0,0,null, subSize);
//
//    System.out.println("end;");
//    out.println("end;");
	String fileName = this.outFile2;
	String binaryFileName [] = fileName.split("[.]");
	File f = new File(binaryFileName[0]+"."+binaryFileName[1]);
	ImageIO.write(img, binaryFileName[1].toUpperCase(), f);

}
private void loadObjFileWithHeights(String heightFile2, double[] fa2, int subSize) throws IOException {
	// TODO Auto-generated method stub
	int rawSize = 512;
	int size = rawSize;
	int size2 = subSize + 1;
	size++;
	double     scale = subSize/rawSize;
	ByteBuffer myByteBuffer = ByteBuffer.allocate(size*size*2);
	myByteBuffer.order(ByteOrder.BIG_ENDIAN);

	ShortBuffer myShortBuffer = myByteBuffer.asShortBuffer();
	
	ByteBuffer myByteBuffer2 = ByteBuffer.allocate(size*size*2);
	myByteBuffer2.order(ByteOrder.BIG_ENDIAN);

	ShortBuffer myShortBuffer2 = myByteBuffer2.asShortBuffer();
	
	for (int i = 0;i <rawSize;i++){
		for (int j = rawSize;j >= 0;j--){
			short s =(short)(64.0 + 127.0 * (fa2[(int) ((i*size2*scale)+j*scale)]));
			System.out.println("s="+s);

			myShortBuffer.put(s);
			myShortBuffer2.put((short)(256-s));
		}
	}

	FileChannel out = new FileOutputStream(heightFile2).getChannel();
	out.write(myByteBuffer);
	out.close();
	FileChannel out2 = new FileOutputStream(heightFile2.split("[.]")[0]+"2"+"."+heightFile2.split("[.]")[1]).getChannel();
	out2.write(myByteBuffer2);
	out2.close();
	System.out.println("Finished");
}


private void loadObjFileWithVertsAndNormals(String string, double[] fa, int subSize) throws IOException {
	// TODO Auto-generated method stub
	int size = subSize;
	size++;
	double     inc = 2.0 / subSize;
	String oText;
	
	double maxx, minx, maxz, minz;
	maxx = minx = maxz = minz =0.0;
	BufferedReader br = new BufferedReader(new FileReader(string));
	try {
	    String line = br.readLine();

	    while (line != null) {
	    	String[] parts = line.split(" ");
	    	if (parts[0].equals("v")){
	    		double x = Double.parseDouble(parts[1]);
	    		double y = Double.parseDouble(parts[2]);
	    		double z = Double.parseDouble(parts[3]);
			    maxx = x + 0.00001 > maxx ? x + 0.00001 : maxx;
			    minx = x - 0.00001  < minx ? x - 0.00001 : minx;    		
			    maxz = z + 0.00001 > maxz ? z + 0.00001 : maxz;
			    minz = z - 0.00001  < minz ? z - 0.00001 : minz;    		
	    	}
	        line = br.readLine();
	    }
	} finally {
	    br.close();
	}	
	
	br = new BufferedReader(new FileReader(string));
	try {
	    StringBuilder sb = new StringBuilder();
	    String line = br.readLine();
	    ArrayList<double[]> alda=new ArrayList<double[]>();
	    int locationOfNormal = 0;
	    while (line != null) {
	    	String[] parts = line.split(" ");
	    	if (parts[0].equals("v")){
	    		double x = Double.parseDouble(parts[1]);
	    		double y = Double.parseDouble(parts[2]);
	    		double z = Double.parseDouble(parts[3]);
	    		double xx = (( 2.0 * x / (maxx - minx)));
	    		double zz = (( 2.0 * z / (maxz - minz)));
	    		int j = (int) ((1.0+xx)*subSize/2.0);
	    		int i = (int) ((1.0+zz)*subSize/2.0);
	    		j = j >= subSize ? subSize -1:j;
	    		i = i >= subSize ? subSize -1:i;
	    		j = j <= 0 ? 0:j;
	    		i = i <= 0 ? 0:i;
	    		
		        System.out.println("xx"+xx+"zz"+zz+"i"+i+"j"+j);

	    		double yy = (( fa[((i+1)*size)+j] * (maxz - minz))/2.0);

			    double normal[] = genNormal (xx, fa[(i*size)+j], zz,
			    		xx, fa[((i+1)*size)+j], zz+inc,
			    		xx+inc, fa[(i*size+j+1)], zz);
			    alda.add(normal);
		        sb.append("v "+x+" "+yy+" "+z);
		        System.out.println("v "+x+" "+yy+" "+z);
		        sb.append(System.lineSeparator());
	    	}else if(parts[0].equals("vn")){
		        sb.append("vn "+alda.get(locationOfNormal)[0]+" "+
        		alda.get(locationOfNormal)[1]+" "+
        		alda.get(locationOfNormal)[2]);
		        System.out.println("*vn "+alda.get(locationOfNormal)[0]+" "+
		        		alda.get(locationOfNormal)[1]+" "+
		        		alda.get(locationOfNormal)[2]);
		        locationOfNormal++;
		        if (locationOfNormal >= alda.size())locationOfNormal = alda.size()-1;
		        sb.append(System.lineSeparator());
	    	}else{
	    		sb.append(line);
		        System.out.println(line);

	    		sb.append(System.lineSeparator());
	    	}
	        line = br.readLine();
	    }
	    oText = sb.toString();

	} finally {
	    br.close();
	}
	
    FileWriter fileWriter =
            new FileWriter(string);

        // Always wrap FileWriter in BufferedWriter.
        BufferedWriter bufferedWriter =
            new BufferedWriter(fileWriter);

        // Note that write() does not automatically
        // append a newline character.
        bufferedWriter.write(oText);
 
        // Always close files.
        bufferedWriter.close();
	
}


private double craterCrossArray(double coeffs[], double r){
	double sum = 0.0;
	for (int i = coeffs.length -1;i >= 0;i--){
		sum += Math.pow(r, (double)i)*coeffs[i];
	}
	return sum;
}
private void divide (double fa[],int w, int ii, int jj, int subSize, PrintWriter out, Dimension d, Graphics g, BufferedImage img){
	int size = subSize;
	size++;
	double     inc = 2.0 / subSize;
	
	double max, min;
	max = min =0.0;
	for (int i = ii;i < ((ii+w+1)<=subSize?ii+w+1:subSize);i++){
		for (int j = jj; j <((jj+w+1)<=subSize?jj+w+1:subSize);j++){
			double x = -1.0 + j*inc;
			double z = -1.0 + i*inc;
			double temp = fa[(i*size)+j];
		    max = temp + 0.00001 > max ? temp + 0.00001 : max;
		    min = temp - 0.00001  < min ? temp - 0.00001 : min;    		
		}
	}
	
	double xx =  -1.0 + jj*inc;
	double zz =  -1.0 + ii*inc;
	double xxw = -1.0 + (jj+w)*inc+inc;
	double zzw = -1.0 + (ii+w)*inc+inc;
	
	System.out.println("bounding "+w+", "+ii+", "+jj+";");
	bounding ( xx,  min,  zz,
		      xxw,  min,  zz,
		      xxw,  max,  zz, 
		      xx,  max,  zz, 
		      xx,  min,  zzw, 
		     out);

	if (w <= 2){
		for (int i = ii;i < ii+w;i++){
			for (int j = jj; j <jj+w;j++){
				double x = -1.0 + j*inc;
				double z = -1.0 + i*inc;
				
			    double normal[] = genNormal (x, fa[(i*size)+j], z,
					    x, fa[((i+1)*size)+j], z+inc,
					    x+inc, fa[(i*size+j+1)], z);
			    draw3DTriangle (x, fa[(i*size)+j], z,
					    x, fa[((i+1)*size)+j], z+inc,
					    x+inc, fa[(i*size+j+1)], z, 
					    out, normal, d, g,j,i,img,subSize);

			    double normal2[] = genNormal (x, fa[((i+1)*size)+j], z+inc,
					    x+inc, fa[((i+1)*size)+j+1], z+inc,
					    x+inc, fa[(i*size+j+1)], z);

			    draw3DTriangle (x, fa[((i+1)*size)+j], z+inc,
					    x+inc, fa[((i+1)*size)+j+1], z+inc,
					    x+inc, fa[(i*size+j+1)], z, 
					    out, normal2, d, g,j,i,img, subSize);

			}
		}
	}else {
		int halfw = w/2;
		divide (fa,halfw, ii, jj, subSize,out,d,g, img);
		divide (fa,halfw, ii+halfw, jj, subSize, out, d, g, img);
		divide (fa,halfw, ii, jj+halfw, subSize, out, d, g, img);
		divide (fa,halfw, ii+halfw, jj+halfw, subSize, out, d, g, img);
	}
	System.out.println("end;");
	if(out!=null)out.println("end;");

}
class RenderCanvas extends Canvas {

		public RenderCanvas(long s, String of, int el){
			seed = s; outFile = of; edgeLength = el;
		}
		public void paint(Graphics g) {
			Dimension d = this.getSize();
			try {
				PrintWriter out = null;
				if (!outFile.equals("null")){
					out = new PrintWriter(new BufferedWriter(new FileWriter(outFile)));
				}
				r = new Random (seed);

				fa = new double [(edgeLength+1) * (edgeLength+1)];
				fill2DFractArray(
									fa,
									edgeLength,
									Long.parseLong(JTFseed.getText()),
									Double.parseDouble(JTFheightScale.getText()),
									Double.parseDouble(JTFh.getText())
									);
				draw2DFractArrayAsTriangles(Long.parseLong(JTFSeed2.getText()),fa, edgeLength, out, d, g);
				if (out != null){
					out.flush();
					out.close();
				}
				
			} catch (Exception e){
				System.out.println(e);
			}
		}

	}

	public static void main(String args[]){

		System.out.println("");
		LandscapeGen lg;

		for(int i = 0; i < args.length; i++){
			System.out.print("args["+i+"] = \""+args[i]+"\", ");
		}
		if(args.length != 8){
			System.out.println("usage: java -classpath bin LandscapeGen seed out_file edge_length type_of_scape(moon,earth) crater_seed normals_image obj_file raw_file");
		}else{
			lg = new LandscapeGen("Landscape Gen", Long.parseLong(args[0]), args[1], Integer.parseInt(args[2]), args[3],Long.parseLong(args[4]), args[5], args[6],args[7]);
		}
	}

}