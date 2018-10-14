import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;

public class my_Bin implements PlugInFilter {
	ImagePlus imp;

	public int setup(String arg, ImagePlus imp) 
        {
		this.imp = imp;
		return DOES_ALL;
	}

	public void run(ImageProcessor ip)
        {
		
		byte[] pixels = (byte[])ip.getPixels();
		int width = ip.getWidth();
		Rectangle r = ip.getRoi();
		int offset,pix, i;
		for (int y=r.y; y<(r.y+r.height); y++)
                {
                     offset = y*width;
                     for (int x=r.x; x<(r.x+r.width); x++)
                     {
                        i = offset + x;
                        //-----  casting a byte variable to a int variable
                        pix = pixels[i] & 0xff ;
                        if(pix > 200)
                           pixels[i] = (byte)(255);
                        else
                           pixels[i] = (byte)(0);
                     }
                }

	}

}
