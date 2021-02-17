#declare CameraBlur = 0; //turns camera blur on and off, 1 = on
#declare AreaLight = 0; //turns softening of shadows on and off, 1 = on

#declare Ni = ~~Ni~~;
#declare Nj = ~~Nj~~;
#declare Nk = ~~Nk~~;
#declare nslice = ~~nslice~~;
//#declare nslice = clock;
#declare Num = str(nslice,0,0);
#warning concat("Num is:", Num,"\n")

#include "colors.inc"

global_settings {
  //max_trace_level 5
  assumed_gamma 2.2
 // radiosity {
//    pretrace_start 0.08
//    pretrace_end   0.01
//    count 35
//    nearest_count 5
//    error_bound 1.8
//    recursion_limit 2
//    low_error_factor .5
   // gray_threshold 0.0
  //  minimum_reuse 0.015
//    brightness 1
//    adc_bailout 0.01/2
  //}
}

#declare lookat = <~~lookat1~~, ~~lookat2~~, nslice>;
#declare focalpoint = lookat;

//focal blur camera
camera {
  location lookat - 880*x
  look_at   lookat
  angle 5
  right     x*image_width/image_height
  #if (CameraBlur)
  	aperture 3          // [0...N] larger is narrower depth of field (blurrier)
  	blur_samples 1000        // number of rays per pixel for sampling
  	focal_point focalpoint    // point that is in focus <X,Y,Z>
  	confidence 0.95           // [0...<1] when to move on while sampling (smaller is less accurate)
  	variance 1/200            // [0...1] how precise to calculate (smaller is more accurate)
  #end
}

//light source
light_source {
  0*x                 // light's position (translated below)
  color rgb 1.0       // light's color
  #if (AreaLight)
  	area_light
  	<0, 0, 200> <200, 0, 0> // lights spread out across this distance (x * z)
  	4, 4                // total number of lights in grid (4x*4z = 16 lights)
  	adaptive 3          // 0,1,2,3...
  	jitter              // adds random softening of light
  	circular            // make the shape of the light circular
  	orient              // orient light
  #end
  translate 200*<-1, .1,  .1>
}


// set a color of the background (sky)
// background { color Black}
background { color srgbt <0.0, 0.0, 0.0, 1.0> }

#declare FOVrad = .3;

#declare FOVplane = union{
//	sphere{<.5,1,.5>,FOVrad}
//	sphere{<16.5,1,.5>,FOVrad}
//	sphere{<16.5,1,16.5>,FOVrad}
//	sphere{<.5,1,16.5>,FOVrad}
	cylinder{<.5,1,.5>,<Ni+.5,1,.5>,FOVrad}
	cylinder{<Ni+.5,1,.5>,<Ni+.5,1,Nj+.5>,FOVrad}
	cylinder{<Ni+.5,1,Nj+.5>,<.5,1,Nj+.5>,FOVrad}
	cylinder{<.5,1,Nj+.5>,<.5,1,.5>,FOVrad}
}

//object{FOVplane pigment{ color Black} no_shadow no_reflection}

#declare F_ODF = finish {
    	ambient .4
    	specular .4
  	}
  	
#declare T_ODF = texture {
 	finish{F_ODF}
    	pigment{rgb <1,1,1>}
  	}

#declare ODF = mesh2 {
  vertex_vectors {
		#declare Fname = concat("ODF_verts_", Num, ".txt");
		#fopen XYZdat Fname read
		#read (XYZdat,DatLength)
		DatLength,
		#declare Count = 0;	
		#while (Count<DatLength-1)
			#read (XYZdat,Xval,Yval,Zval)
			<Xval, Zval, Yval>,
			#declare Count = Count+1;	
		#end
		#read (XYZdat,Xval,Yval,Zval)
		<Xval, Zval, Yval>
		#fclose XYZdat
  }
  
  	normal_vectors {
		#declare Fname = concat("ODF_norms_",Num,".txt");
		#fopen XYZdat Fname read
		#read (XYZdat,DatLength)
		DatLength,
		#declare Count = 0;	
		#while (Count<DatLength-1)
		#read (XYZdat,Xval,Yval,Zval)
			<Xval, Zval, Yval>,
			#declare Count = Count+1;	
		#end
		#read (XYZdat,Xval,Yval,Zval)
		<Xval, Zval, Yval>
		#fclose XYZdat
	}
	  
	  texture_list {
		#declare Fname = concat("ODF_c~~char_r~~_",Num,".txt");
		#fopen XYZdat Fname read
		#read (XYZdat,DatLength)
		DatLength,
		#declare Count = 0;	
		#while (Count<DatLength-1)
		#read (XYZdat,Xval,Yval,Zval)
			texture{finish{F_ODF} pigment{rgb <Xval, Yval, Zval>}}
			#declare Count = Count+1;	
		#end
		#read (XYZdat,Xval,Yval,Zval)
		texture{finish{F_ODF} pigment{rgb <Xval, Yval, Zval>}}
		#fclose XYZdat
	   }

	face_indices {
		#declare Fname = concat( "ODF_tri_",Num,".txt");
		#fopen XYZdat Fname read
		#read (XYZdat,DatLength)
		DatLength,	
		#declare Count = 0;	
		#while (Count<DatLength-1)
			#read (XYZdat,Xval,Yval,Zval)
			<Xval, Yval, Zval>, Xval, Yval, Zval
			#declare Count = Count+1;
		#end
		#read (XYZdat,Xval,Yval,Zval)
		<Xval, Yval, Zval>, Xval, Yval, Zval
		#fclose XYZdat
	}  
}

object{ODF scale <1,1,1> rotate 0*z rotate -0*y no_shadow no_reflection}



