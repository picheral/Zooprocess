import ij.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.gui.*;
import ij.measure.*;
import ij.text.*;
import java.util.*;

/*
Convex Hull Plus by Gabriel Landini, 12/Sep/2004.  G.Landini at bham. ac. uk

Minimum bounding circle as suggested by:
From:         Xavier Draye <draye@ECOP.UCL.AC.BE>
Sender:       ImageJ Interest Group <IMAGEJ@LIST.NIH.GOV>
Subject:      Re: Bounding circle

If you are still interested in getting the smallest circle enclosing all
your points, I think the following would do the job:

REASONING: the number of points that will be exactly on the circle border
will be at least 2.
-If there are 2 points on the circle, the circle center will be the
midpoint of the segment joining the two points.
-If there are 3 points on the circle, we can find for each side of the
triangle the unique line that contains the midpoint the side segment and is
orthogonal to the segment. There will be 3 such lines, and the center of
the circle is the intersection of any two of them.
-If there are more than three points, then the situation reduces to the 3
points (taking any three of them, although there will be choices of three
that may lead to greater numerical errors). Actually, I am not sure the >3
points case will not be already solved by the 2 or 3 points case.

ONE implementation of this would be:
Let us calculate the distance between all possible pairs of points. We then
take the pair of points (A & B) that are the most distant. Let P be the
midpoint of segment AB. Then find the point C which is neither A nor B and
whose distance to P is the greatest.
-If the distance CP is smaller than the distance PA (or PB), we are in the
2 points case and P is the center of the circle we are looking for (I will
call it the 2-points case circle).

-If the distance CP is greater than the distance PA (or PB), we are in
the >2 points case. 

------
Note: Adrian Daerr pointed out that the original suggestion (below) does 
not always work to find the minimal bounding circle:
------

  Geometry tells us that A and B must also be on the (new) circle border 
  (because A and B are the most distant points). Also, we know that the 
  third point we are looking for is outside the circle of the 2 points 
  case. It would seem logical to me that C is that third point which
  we can check by calculating the new circle and checking all points are
  within it.


So, instead of using the test above, one can search for the smallest 
osculating circle defined by 3 points in the convex hull that contains 
all convex hull points (and therefore all points in the set since these
are contained within the convex hull).

*/

// 26/12/2004 returns floating values.
// 22/5/2005 changed test for the case  of 3 points defining the osculating circle.

public class Convex_Hull_esmer implements PlugInFilter, Measurements {
	ImagePlus imp;
	int counter=0;
	public int setup(String arg, ImagePlus imp) {
		ImageStatistics stats;
		stats=imp.getStatistics();
		if (stats.histogram[0]+stats.histogram[255]!=stats.pixelCount){
			IJ.error("8-bit binary image (0 and 255) required.");
			return DONE;
		}

	
		if (arg.equals("about"))
			{showAbout(); return DONE;}
		this.imp = imp;
		return DOES_8G;
	}

	public void run(ImageProcessor ip) {
		GenericDialog gd = new GenericDialog("Convex Hull Plus 1.2", IJ.getInstance());
		//gd.addMessage("Set must be white!");
		String [] mode={"Convex Hull selection", "Minimal Bounding Circle selection", "Draw Convex Hull" ,"Draw Minimal Bounding Circle","Draw both"};
		gd.addChoice("Mode", mode, mode[0]);
		gd.addCheckbox("White particles",true);
		gd.addCheckbox("Log results",true);

		gd.showDialog();
		if (gd.wasCanceled())
			return;
		String myMode = gd.getNextChoice();
		boolean white = gd.getNextBoolean();
		boolean logr = gd.getNextBoolean();

		int i, j, k=0, m, colour=0;

		if (white) 
			colour=255;
		//IJ.run("Set Scale...", "distance=0 known=1 pixel=1 unit=pixels"); // OK?
		imp.setCalibration(null); // OK?

		for (j=0;j<imp.getHeight();j++){
			for (i=0;i<imp.getWidth();i++){
				if (ip.getPixel(i,j)== colour)
					counter++;
			}
		}

		int[] x = new int[counter+1];
		int[] y = new int[counter+1];

		for (j=0;j<imp.getHeight();j++){
			for (i=0;i<imp.getWidth();i++){
				if (ip.getPixel(i,j) == colour){
					x[k] = i;
					y[k] = j;
					k++;
				}
			}
		}

		//cnvx hull
		int n=counter, min = 0, ney=0, px, py, h, h2, dx, dy, temp, ax, ay;
		double minangle, th, t, v, zxmi=0;

		for (i=1;i<n;i++){
			if (y[i] < y[min]) 
				min = i;
		}

		temp = x[0]; x[0] = x[min]; x[min] = temp;
		temp = y[0]; y[0] = y[min]; y[min] = temp;
		min = 0;

		for (i=1;i<n;i++){
			if (y[i] == y[0]){
				ney ++;
				if (x[i] < x[min]) min = i;
			}
		}
		temp = x[0]; x[0] = x[min]; x[min] = temp;
		temp = y[0]; y[0] = y[min]; y[min] = temp;
		ip.setColor(127);

		//first point x(0), y(0)
		px = x[0];
		py = y[0];

		min = 0;
		m = -1; 
		x[n] = x[min];
		y[n] = y[min];
		if (ney > 0) 
			minangle = -1;
		else
			minangle = 0;

		while (min != n+0 ){
			m = m + 1;
			temp = x[m]; x[m] = x[min]; x[min] = temp;
			temp = y[m]; y[m] = y[min]; y[min] = temp;

			min = n ; //+1
			v = minangle;
			minangle = 360.0;
			h2 = 0;

			for (i = m + 1;i<n+1;i++){
				dx = x[i] - x[m];
				ax = Math.abs(dx);
				dy = y[i] - y[m];
				ay = Math.abs(dy);
  
				if (dx == 0 && dy == 0) 
					t = 0.0;
				else 
					t = (double)dy / (double)(ax + ay);
  
				if (dx < 0)
					t = 2.0 - t;
				else {
					if (dy < 0)
						t = 4.0 + t;
				}
				th = t * 90.0;

				if(th > v){
					if(th < minangle){
						min = i;
						minangle = th;
						h2 = dx * dx + dy * dy;
					}
					else{
						if (th == minangle){
							h = dx * dx + dy * dy;
							if (h > h2){
								min = i;
								h2 = h;
							}
						}
					}
				}
			}
			if (myMode.equals("Draw Convex Hull") || myMode.equals("Draw both"))
				ip.drawLine(px,py,x[min],y[min]);	
			px = x[min];
			py = y[min];
			zxmi = zxmi + Math.sqrt(h2);
		}
		m++;
			
		int[] hx = new int[m];// ROI polygon array 
		int[] hy = new int[m];

		for (i=0;i<m;i++){
			hx[i] =  x[i]; // copy to new polygon array
			hy[i] =  y[i];
		}
		
		if (myMode.equals("Convex Hull selection"))
			imp.setRoi(new PolygonRoi(hx, hy, hx.length, 2)); // roi.POLYGON

		if (logr){
			IJ.log("Hull_points: "+ m);
			IJ.log("Hull_length: "+(float)zxmi);
		}

		// get the edges between points
		double [] d = new double [(m *(m-1))/2]; // edge lenght
		int[] p1 = new int [(m *(m-1))/2]; // point 1
		int[] p2 = new int [(m *(m-1))/2]; // point 2

		k=0;
		for (i=0;i<m-1;i++){
			for (j=i+1;j<m;j++){
				d[k]= Math.sqrt(Math.pow(hx[i]-hx[j], 2.0) + Math.pow(hy[i]-hy[j], 2.0));
				p1[k]=i;
				p2[k]=j;
				k++;
			}
		}
		k--;

		// sort the distances
		boolean sw = true;
		double tempd;
		int centre, cw, pc, p3;
		double tt, tttemp, radius, cx, cy;
		
		while (sw){
			sw=false;
			for(i=0;i<k-1;i++){
				if (d[i]<d[i+1]){
					tempd = d[i]; d[i] = d[i+1]; d[i+1] = tempd;
					temp = p1[i]; p1[i] = p1[i+1]; p1[i+1] = temp;
					temp = p2[i]; p2[i] = p2[i+1]; p2[i+1] = temp;
					sw = true;
				}
			}
		}

		//IJ.log("1:"+d[0]+" "+p1[0]+" "+p2[0]);
		//IJ.log("2:"+d[1]+" "+p1[1]+" "+p2[1]);
		radius=d[0]/2.0;

		cx=(hx[p1[0]]+hx[p2[0]])/2.0;
		cy=(hy[p1[0]]+hy[p2[0]])/2.0;

		// find largest distance from point c
		//sw=false;
		p3=-1;
		tt=radius;
		for (i=0;i<m;i++){
			tttemp=Math.sqrt(Math.pow(hx[i]-cx, 2.0) + Math.pow(hy[i]-cy, 2.0));
			if(tttemp>tt){
				tt=tttemp;
				//IJ.log("Largest from c: "+Math.sqrt(Math.pow(hx[i]-cx, 2.0) + Math.pow(hy[i]-cy, 2.0)));
				p3=i;
			}
		}

		if (p3>-1){
			// 3 or more- point circle
			//IJ.log("p3:"+p3);
			//IJ.log("osculating circle p1, p2, p3");
			double [] op1 = new double [2];
			double [] op2 = new double [2];
			double [] op3 = new double [2];
			double [] circ = new double [3];
			double tD=Double.MAX_VALUE;
			int tp1=0, tp2=0, tp3=0, z;
			
			// GL new test: find the smallest osculating circle that contains all convex hull points
			for (i=0; i<m-2; i++){
				for (j=i+1; j<m-1; j++){
					for (k=j+1; k<m; k++){
						op1[0]=hx[i];
						op1[1]=hy[i];
						op2[0]=hx[j];
						op2[1]=hy[j];
						op3[0]=hx[k];
						op3[1]=hy[k];
						osculating(op1, op2, op3, circ);
						// IJ.log(""+i+ " "+j+" "+k+"   "+circ[2]);
						// store a large dummy radius
						if (circ[2]>0){
							sw=true;
							for (z=0;z<m;z++){
								tttemp=(float)Math.sqrt(Math.pow(hx[z]-circ[0], 2.0) + Math.pow(hy[z]-circ[1], 2.0));
								if(tttemp>circ[2]){
									sw=false; // points are outside it
									break; // don't check any more points
								}
							}
							if(sw){ //no CH points outside
								// store radius & coordinates
								// IJ.log(""+i+ " "+j+" "+k+"   "+circ[2]);
								if (circ[2]<tD){
									tp1=i;
									tp2=j;
									tp3=k;
									tD=circ[2];
								}
							}
						}
					}
				}
			}
			op1[0]=hx[tp1];
			op1[1]=hy[tp1];
			op2[0]=hx[tp2];
			op2[1]=hy[tp2];
			op3[0]=hx[tp3];
			op3[1]=hy[tp3];
			//IJ.log("Solved for:"+tp1+ " "+tp2+" "+tp3);
			osculating(op1, op2, op3, circ);
			radius=circ[2];
			if (logr)
				IJ.log("Bounding circle (3 points) x: "+(float) circ[0]+", y: "+(float) circ[1]+", radius: "+ radius);
			if (myMode.equals("Minimal Bounding Circle selection")){
				if (circ[2]>0)
					IJ.makeOval((int) Math.floor((circ[0]-circ[2])+.5), (int) Math.floor((circ[1]-circ[2])+.5),(int)Math.floor((radius*2)+.5), (int)Math.floor((radius*2)+.5));	
			}
			if (myMode.equals("Draw Minimal Bounding Circle") || myMode.equals("Draw both")){
				if (circ[2]>0){
					IJ.makeOval((int) Math.floor((circ[0]-circ[2])+.5), (int) Math.floor((circ[1]-circ[2])+.5), (int)Math.floor((radius*2)+.5), (int)Math.floor((radius*2)+.5));	
					IJ.run("Unlock Image");
					IJ.run("Draw");
				}
			}
		}
		else{
			//2-point circle centred at cx, cy radius
			if (myMode.equals("Minimal Bounding Circle selection"))
				IJ.makeOval((int)Math.floor(cx-radius+.5), (int)Math.floor(cy-radius+.5), (int)(Math.floor(d[0]+.5)),(int)(Math.floor(d[0]+.5)));
			if (logr)
				IJ.log("Bounding circle (2 points) x: "+(float) cx+", y: "+(float) cy+", radius: "+(float) radius);
			if (myMode.equals("Draw Minimal Bounding Circle") || myMode.equals("Draw both")){
				IJ.makeOval((int)(cx-radius), (int)(cy-radius), (int)(Math.floor(d[0]+.5)),(int)(Math.floor(d[0]+.5)));
				IJ.run("Unlock Image");
				IJ.run("Draw");
			}
		}
	}

	
	void osculating( double [] pa, double [] pb, double [] pc, double [] centrad){
		// returns 3 double values: the centre (x,y) coordinates & radius
		// of the circle passing through 3 points pa, pb and pc
		double a, b, c, d, e, f, g;

		if ((pa[0]==pb[0] && pb[0]==pc[0]) || (pa[1]==pb[1] && pb[1]==pc[1])){ //colinear coordinates
			centrad[0]=0; //x
			centrad[1]=0; //y
			centrad[2]=-1; //radius
			return;
		}

		a = pb[0] - pa[0];
		b = pb[1] - pa[1];
		c = pc[0] - pa[0];
		d = pc[1] - pa[1];

		e = a*(pa[0] + pb[0]) + b*(pa[1] + pb[1]);
		f = c*(pa[0] + pc[0]) + d*(pa[1] + pc[1]);

		g = 2.0*(a*(pc[1] - pb[1])-b*(pc[0] - pb[0]));
		//  If g is 0 then the three points are colinear and no finite-radius
		//  circle through them exists. Return -1 for the radius. Somehow this does not
		// work as it should (representation of double number?), so it is trapped earlier..

		if (g==0.0){
			centrad[0]=0; //x
			centrad[1]=0; //y
			centrad[2]=-1; //radius
		}
		else {//return centre and radius of the circle
			centrad[0] = (d * e - b * f) / g;
			centrad[1] = (a * f - c * e) / g;
			centrad[2] = (float)Math.sqrt(Math.pow((pa[0] - centrad[0]),2) + Math.pow((pa[1] - centrad[1]),2));
		}
	}
	
	
	
	void showAbout() {
		IJ.showMessage("About Convex Hull Plus...",
		"Convex Hull Plus by Gabriel Landini,  G.Landini@bham.ac.uk\n"+
		"ImageJ plugin for calculating the Convex Hull and\n"+
		"the Minimal Bounding Circle of a point set.");
	}
}



