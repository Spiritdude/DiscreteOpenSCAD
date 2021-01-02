// == DiscreteOpenSCAD, An Extension of OpenSCAD ==
//    written by Rene K. Mueller <spiritdude@gmail.com>
//
// License: MIT
//
// Description:
//    Reimplements some of the basic solids in order to maintain
//       accessibility of the mesh before it's transformed into non-accessible objects
//    This approach allows low-level manipulation of vertices, e.g. extruding from one polygon to another
//    See also https://xyzdims.com/discrete-openscad/
//
// History:
// 2021/01/02: 0.0.5: release code finally https://github.com/Spiritdude/DiscreteOpenSCAD
// 2018/12/26: 0.0.4: 3d: dm_extrude_rotate(p,n,start=0,end=360) which allows simple dm_torus()
// 2018/12/25: 0.0.3: 2d: added dp_translate(), dp_scale(), dp_rotate(), dp_bounds() and dp_center(); 3d: added dm_sphere(), dm_cylinder() and dm_cube(), and dm_translate(), dm_scale(), dm_rotate(), dm_bounds() and dm_center() as well
// 2018/12/24: 0.0.2: dp_morph() as part of dm_extrude() uses dp_nearest() to find nearby point of two polygons to morph smoothly, morph examples added
// 2018/12/23: 0.0.1: Point Functions, Polygon Functions and Mesh Functions added, code not yet released


_aoff = 0;     // dpf_npolygon() & dpf_ellipse() is used in a lot low-level polygon creating functions where
               //    _aoff is the angle offset
               //    don't change unless you know what you are doing: it changes the starting point of the polygons
               //       _aoff = 0, top(y=1)/center(x=0) first edge of any npolygon created polygon
               
// ---------------------------------------------------------------------------------------------------------------------
// == npolygon 
function dpf_npolygon(r=1,s=3,p,a=0,aoff=0,poff=0,start=0,end=360) = 
   /*
      Example:
      s = 4          0-90, 90-180, 180-270, 270-360
      p = 0.00       ns=0 f=0 p1=0  p2=90
          0.24    -> ns=0 f=1 p1=0  p2=90
          0.25    -> ns=1 f=0 p1=90 p2=180
          0.50    -> ns=2 f=0 p1=180 p2=270
          0.75    -> ns=3 f=0 p1=270 p2=360
    */
   let(
      ac = (end-start)+start,
      ns = floor((p+poff)%1*s), 
      p1 = (ac/s)*ns+aoff, 
      p2 = (ac/s)*(ns+1)+aoff, 
      f = ((p+poff)%1*s*s)%s/s            // fine-grained
   )   
   [ 
      r*lerp(sin(p1),sin(p2),f),r*lerp(cos(p1),cos(p2),f), 
   ] * [[cos(a+_aoff),-sin(a+_aoff)],[sin(a+_aoff),cos(a+_aoff)]];

function dp_npolygon(r,s,n,a=0,aoff=0,poff=0,start=0,end=360) = [
   let(n=n>0?n:s*8)
   for(i=[0:n-1]) 
      dpf_npolygon(r,s,i/n,a=a,aoff=aoff,poff=poff,start=start,end=end)
];

module dd_npolygon(r,s,n=12,a=0) {
   polygon(dp_npolygon(r,s,n=n,a=a));
}

// ---------------------------------------------------------------------------------------------------------------------
// == ellipse
function dpf_ellipse(r1=2,r2=1,p,a=0) = [ r1*sin(360*p),r2*cos(360*p) ] * 
      [[cos(a+_aoff),-sin(a+_aoff)],[sin(a+_aoff),cos(a+_aoff)]];

function dp_ellipse(r1=2,r2=1,n=32,a=0) = [
      for(i=[0:n-1]) 
         dpf_ellipse(r1,r2,i/n,a=a)
   ];

module dd_ellipse(r1=2,r2=1,n=12,a=0) {
   polygon(dp_ellipse(d1>0?d1/2:r1,d2>0?d2/2:r2,n=n,a=a));
}

// == arc
function dp_arc(r,start=0,end=360,n=32,a=0,type="pie") = [
      for(i=[start:(end-start)/n:end]) 
         dpf_ellipse(r,r,i/360,a=a)
      ,
      type=="pie"?[0,0]:[0,0]
   ];

module dd_arc(r,start=0,end=360,n=12,a=0) {
   polygon(dp_arc(d>0?d/2:r,start,end,n=n,a=a));
}

// ---------------------------------------------------------------------------------------------------------------------
// == circle
function dpf_circle(r=1,a=0) = dpf_ellipse(d>0?d/2:r,d>0?d/2:r,a=a);

function dp_circle(r=1,n=32,a=0) = dp_ellipse(d>0?d/2:r,d>0?d/2:r,n=n,a=a);

module dd_circle(r=1,n=12,a=0) {
   dd_ellipse(d>0?d/2:r,d>0?d/2:r,n=n,a=a);
}

// ---------------------------------------------------------------------------------------------------------------------
// == rectangle
function dpf_rectangle(w=2,h=1,p,a=0) = let(p=dpf_npolygon(sqrt(2)/2,4,p,a=a,aoff=-45,poff=1/8)) [p[0]*w,p[1]*h];
function dpf_rectangle2(w=2,h=1,p,a=0,poff=0) =
   let(p0=(p+poff)%1)
   (p0 < 0.25 ? 
      [ lerp(-w/2,w/2,p0*4),-h/2 ] :
      p0 < 0.5 ? 
         [ w/2, lerp(-h/2,h/2,(p0-0.25)*4) ] :
            p0 < 0.75 ? 
               [ lerp(w/2,-w/2,(p0-0.5)*4), h/2 ] :
               [ -w/2, lerp(h/2,-h/2,(p0-0.75)*4) ]) *
   [[cos(a),-sin(a)],[sin(a),cos(a)]];
            
function dp_rectangle(w=2,h=1,n=8,a=0) = [ 
   let(n=n>3?n:12)
   for(i=[0:n-1]) 
      dpf_rectangle(w,h,i/n,a=a) 
];

//function dp_rectangle(w=2,h=1) = [
//   [-w/2,-h/2],
//   [w/2,-h/2],
//   [w/2,h/2],
//   [-w/2,h/2]
//];

module dd_rectangle(w=2,h=1,a=0) {
   polygon(dp_rectangle(w,h));
}

// ---------------------------------------------------------------------------------------------------------------------
// == square
function dpf_square(w,f,a=0) = dpf_rectangle(w,w,f,a=a);

function dp_square(w=1,n=8,a=0) = dp_rectangle(w,w,n=n,a=a);

module dd_square(w=1,n=12,a=0) {
   polygon(dp_square(w,n=n,a=a));
}

// ---------------------------------------------------------------------------------------------------------------------
// == polygon manipulation

function dp_subdivide(p,n=2) = [
   let(m=len(p))
   for(i=[1:m])
      for(j=[0:n-1]) [
         for(k=[0:1]) 
            lerp(p[i-1][k],p[i%m][k],j/n)
      ]
];

function dp_distort(p,s=1) = [
   let(m=len(p))
   for(i=[0:m-1]) [
      for(k=[0:1]) 
         p[i][k]+(1-2*rand())*s
   ]   
];

function dp_morph(a,b,p) = let(n=len(a),m=len(b),q=n>m?n:m) [
   for(i=[0:q-1]) 
      //let(near=dp_nearest(b,a[i])[0])
      [
         //lerp(a[i].x,b[near].x,p),
         //lerp(a[i].y,b[near].y,p)
         lerp(a[i*n/q].x,b[i*m/q].x,p),
         lerp(a[i*n/q].y,b[i*m/q].y,p)
      ]
];

function dp_nearest(a,p,i=0) = let(     // experimental
      d = sqrt( (a[i].x-p.x)*(a[i].x-p.x) + (a[i].y-p.y)*(a[i].y-p.y) ),
      d1 = [i,d],
      d2 = i < len(a)-1 ? dp_nearest(a,p,i+1) : [i,d]
   ) d1[1] < d2[1] ? d1 : d2;

function dp_scale(p,s) = [
   let(s=len(s)>0?s:[s,s])
   for(i=[0:len(p)-1]) [
      for(j=[0:1]) 
         p[i][j] * s[j]
   ]
];

function dp_rotate(p,a=0) = [
   let(mr=[[cos(a),-sin(a)],[sin(a),cos(a)]])
   for(i=[0:len(p)-1])
      p[i] * mr
];

function dp_translate(p,t=0) = [
   let(t0=len(t)>0?t:[t,t])
   for(i=[0:len(p)-1]) 
      p[i] + t0
];

function dp_bounds(p,i=0) =
   let(
      a = [[p[i].x,p[i].x],[p[i].y,p[i].y]],
      b = i < len(p)-1 ? dp_bounds(p,i+1) : a,
      minx = a.x[0] < b.x[0] ? a.x[0] : b.x[0],
      maxx = a.x[1] > b.x[1] ? a.x[1] : b.x[1],
      miny = a.y[0] < b.y[0] ? a.y[0] : b.y[0],
      maxy = a.y[1] > b.y[1] ? a.y[1] : b.y[1]
   )
   [ [minx,maxx],[miny,maxy] ];

function dp_center(p,center=[true,true]) = 
   let(b=dp_bounds(p))
   dp_translate(p,[center[0]?(-b.x[0]-b.x[1])/2:0,center[1]?(-b.y[0]-b.y[1])/2:0]);

// ---------------------------------------------------------------------------------------------------------------------
// == extrude

function dm_extrude(a,b,n=12,h=1) = 
   let(
      h = h > 0 ? h : 1,
      //ab = len(a) > len(b) ? [a,b] : [b,a],
      //a = ab[0],
      //b = ab[1],
      al = len(a),
      bl = len(b),
      m = al > bl ? al : bl,
      //m = al, //len(a),
      off = n*m,
      nmb = [                        // for computational reasons, create near map
         for(i=[0:len(a)-1])
            dp_nearest(b,a[i])[0]
      ],
      nma = [                        // for computational reasons, create near map
         for(i=[0:len(b)-1])
            dp_nearest(a,b[i])[0]
      ],
      // points
      p = [ 
         for(k=[0:n]) 
            //for(i=[0:k<n/2?len(a)-1:len(b)-1])
            //for(i=[0:len(a)-1])
            for(i=[0:m-1])
               //let(near=dp_nearest(b,a[i])[0])
               //[ lerp(a[i].x,b[i].x,k/n), lerp(a[i].y,b[i].y,k/n), h*(k/n) ],
               //[ lerp(a[i].x,b[near].x,k/n), lerp(a[i].y,b[near].y,k/n), h*(k/n) ],
                  //k < n/2 ? 
                  //   [ lerp(a[i].x,b[nmb[i]].x,k/n), lerp(a[i].y,b[nmb[i]].y,k/n), h*(k/n) ],
                  //:
                  //   [ lerp(a[nma[i]].x,b[i].x,k/n), lerp(a[nma[i]].y,b[i].y,k/n), h*(k/n) ],
                  [ lerp(a[i*al/m].x,b[i*bl/m].x,k/n), lerp(a[i*al/m].y,b[i*bl/m].y,k/n), h*(k/n) ],
                  
         [0,0,0],
         [0,0,h]
      ],
      // faces
      f0 = [ 
         for(i=[0:m-1]) 
            [
               i,
               (i+1)%m,
               len(p)-2
            ]
      ],
      f1 = [ 
         for(k=[1:n]) 
            for(i=[1:m]) 
               [
                  k*m+i-1,
                  k*m+(i%m),
                  (k-1)*m+i%m
               ]
      ],
      f2 = [ 
         for(k=[1:n]) 
            for(i=[1:m]) 
               [
                  k*m+i-1,
                  (k-1)*m+i%m,
                  (k-1)*m+i-1
               ]
      ],
      f3 = [
         for(i=[0:m-1]) 
            [
               off+i%m,
               len(p)-1,
               off+(i+1)%m,
            ]
      ]
   )
   [
      p,concat(f0,f1,f2,f3)
   ];

function dm_extrude_rotate(p,n=4,a=0,start=0,end=360) = [
   let(m=len(p))
   [  // points
      for(j=[0:n]) 
         let(a=(end-start)/n*j+start+a)
         for(i=[0:m-1]) 
            let(p=p[i],q=[p.x,0]*[[cos(a),-sin(a)],[sin(a),cos(a)]])
            [ q.x, q.y, p.y ]
   ],
   let(m=len(p),n=(end!=360||start!=0)?n+1:n,off=(end!=360||start!=0)?1:0)
   [  // faces
      for(j=[0:n]) 
         for(i=[0:m-1]) 
           [ 
              j*m+i, 
              j*m+(i+1)%m, 
              ((j+1)%(n+off))*m+(i+1)%m 
           ],
      for(j=[0:n]) 
         for(i=[0:m-1]) 
           [ 
              j*m+i, 
              ((j+1)%(n+off))*m+(i+1)%m,
              ((j+1)%(n+off))*m+i
           ]
   ],
];

function dm_extrude_along(a,b,path,n=4) = 
   let(
      path = cm_smooth3(path,n),
      m = len(path)
   )
   dm_merge([
      for(i=[1:m-1])     // each smoothed path point
         let(
            p0 = path[i-1],
            p1 = path[i],
            d = norm(p1,p0)
         )
         dm_translate(p0,
            dm_rotate(
               dm_extrude(a,b,i/(m-1)),
               look_at(p0,p1)
            )
         )
   ]);

function dm_distort(m,s=1,f=1) = [
   // points
   [
      let(l=len(m[0]))
      for(i=[0:l-1]) [
         for(k=[0:2]) 
            //m[0][i][k] + sin(m[0][i][0]*5+1) * cos(m[0][i][1]*5+2) * sin(m[0][i][2]*5 +2) * (1-2*rand()*s)
            //m[0][i][k] + (1-2*rand()*s)
            m[0][i][k] + (1-2*noise3(m[0][i],f))*s
      ]
   ],
   // faces
   m[1]
];

function dm_scale(m,s) = [
   [
      let(s0=len(s)>0?s:[s,s,s])
      for(i=[0:len(m[0])-1]) [
         for(j=[0:2]) 
            m[0][i][j] * s0[j]
      ]
   ],
   m[1]
];

function dm_rotate(m,r) = [
   [
      let(
         r = len(r)>0 ? r : [r,r,r],
         _mr = [
            // reference http://www.songho.ca/opengl/gl_anglestoaxes.html  Rx,Ry,Rz
            //[cos(r.y)*cos(r.z),sin(r.x)*sin(r.y)*cos(r.z)+cos(r.x)*sin(r.z),-cos(r.x)*sin(r.y)*cos(r.z)+sin(r.x)*sin(r.z)],
            //[-cos(r.y)*sin(r.z),-sin(r.x)*sin(r.y)*sin(r.z)+cos(r.x)*cos(r.z),cos(r.x)*sin(r.y)*sin(r.z)+sin(r.x)*cos(r.z)],
            //[sin(r.y),-sin(r.x)*cos(r.y),cos(r.x)*cos(r.y)]
            [cos(r.y)*cos(r.z),-cos(r.y)*sin(r.z),sin(r.y)],
            [sin(r.x)*sin(r.y)*cos(r.z)+cos(r.x)*sin(r.z),-sin(r.x)*sin(r.y)*sin(r.z)+cos(r.x)*cos(r.z),-sin(r.x)*cos(r.y)],
            [-cos(r.x)*sin(r.y)*cos(r.z)+sin(r.x)*sin(r.z),cos(r.x)*sin(r.y)*sin(r.z)+sin(r.x)*cos(r.z),cos(r.x)*cos(r.y)]
         ],
         mr = mt_rotate(r)
      )
      for(i=[0:len(m[0])-1])
         //m[0][i] * _mr
         mt_transform(vec3(m[0][i]),mr)
   ],
   m[1]
];

function dm_translate(m,t) = [
   [
      let(t0=len(t)>0?t:[t,t,t])
      for(i=[0:len(m[0])-1]) [
         for(j=[0:2]) 
            m[0][i][j] + t0[j]
      ]
   ],
   m[1]
];

// ---------------------------------------------------------------------------------------------------------------------
// simple solids

function dm_cube(s=1,center=false) = dm_extrude(
   dp_rectangle(len(s)?s[0]:s,len(s)?s[1]:s,n=8),
   dp_rectangle(len(s)?s[0]:s,len(s)?s[1]:s,n=8),
   h=len(s)?s[2]:s,n=1);

function dm_cylinder(r=1,h=4,n=12,center=false,d=0,r1=undef,r2=undef,d1=undef,d2=undef) = dm_extrude(
   dp_circle(d>0?d/2:(r1!=undef?r1:(d1!=undef?d1/2:r)),n=n),
   dp_circle(d>0?d/2:(r2!=undef?r2:(d2!=undef?d2/2:r)),n=n),
   h=h);

function dm_sphere(r=1,center=false,n=12,d=0) = [
   [ // points
      let(r=d>0?d/2:r)
      for(j=[0:n]) 
         for(i=[0:n-1]) 
            [ sin(360/n*i)*r * sin(180/n*j), cos(360/n*i)*r * sin(180/n*j), cos(180/n*j)*r ]
   ],
   [  // faces
      for(j=[0:n-2]) 
         for(i=[0:n-1]) 
            [ 
               j*n+i, 
               (j+1)*n+i,
               (j+1)*n+(i+1)%n,
            ],
      for(j=[0:n-1]) 
         for(i=[0:n-1]) 
            [ 
               j*n+i, 
               j*n+(i+1)%n,
               (j+1)*n+(i+1)%n,
            ]
   ]
];

function dm_torus(ri=1,ro=4,ni=48,no=48,ai=0,ao=0) = 
   dm_extrude_rotate(dp_translate(dp_npolygon(ri,ni,a=ai),[ro,0]),no,a=ao);

// ---------------------------------------------------------------------------------------------------------------------
// == mesh manipulation

function dm_matrix_mul(m,ma) = [
   [
      for(i=[0:len(m[0])-1]) 
         m[0][i] * ma
   ],
   m[1]
];

function look_at(v0,v1) = 
   let(v=v1-v0)
   [
      0,
      acos(v[2]/norm(v0-v1)),
      atan2(v[1],v[0])
   ];
   
function dm_line(p0,p1,r=1) =
   let(r=d>0?d/2:r,v=p0-p1)
   dm_translate(
      dm_rotate(
         dm_cylinder(r,norm(v)),
         look_at(p0,p1)
      ),
      p0
   );

function min1(a,b) = a < b ? a : b;
function max1(a,b) = a > b ? a : b;
   
function dm_unique_edges(m) = [
            for(i=[0:len(m[1])-1])     // walk through all faces
               let(fn=len(m[1][i]))
               for(j=[0:fn-1]) 
                  //min1(m[1][i][j],m[1][i][(j+1)%fn]),
                  //max1(m[1][i][j],m[1][i][(j+1)%fn])
                  if(m[1][i][j]<m[1][i][(j+1)%fn]) 
                     [m[1][i][j],m[1][i][(j+1)%fn]]
];

function _face_offset(ms,i) = (i < 1  ? 0 : len(ms[i-1][0])) + (i > 0 ? _face_offset(ms,i-1) : 0);

function dm_merge(ms) = [
   [ // points
      // push all points of all meshes
      for(m=[0:len(ms)-1]) 
         for(i=[0:len(ms[m][0])-1])
            ms[m][0][i]
   ],
   let(
      off = [
         for(i=[0:len(ms)-1]) 
            _face_offset(ms,i)
     ]
   )
   [ // faces
      // push all faces with the respective offset per mesh
      for(m=[0:len(ms)-1]) 
         for(i=[0:len(ms[m][1])-1]) [     // all faces
            for(j=[0:2]) 
               ms[m][1][i][j]+off[m]
         ]
   ]
];
            
function dm_wireframe(m,r=1) = 
   // issues:
   //    - when walking through "faces", many edges are shared between faces
   //      so walking faces -> edges will create many double used edges
   //      => avoid duplications
   //      faces = [[0,1,2],[2,3,0],[1,2,4]]
   //              [0,1,2] => [0,1],[1,2],
   //              [2,3,0] => [2,3],[0,3],
   //              [1,2,4] => [1,2]*,[2,4]      *) double
   let(
      edges = dm_unique_edges(m)
   )
   dm_merge([  
      for(i=[0:len(edges)-1]) 
         dm_line(m[0][edges[i][0]],m[0][edges[i][1]],r=r) 
   ]);
   //dm_line(m[0][edges[0][0]],m[0][edges[0][1]],r=r);

// ---------------------------------------------------------------------------------------------------------------------

function dm_bounds(m,i=0) =
   let(
      a = [[m[0][i].x,m[0][i].x],[m[0][i].y,m[0][i].y],[m[0][i].z,m[0][i].z]],
      b = i < len(m[0])-1 ? dm_bounds(m,i+1) : a,
      minx = a.x[0] < b.x[0] ? a.x[0] : b.x[0],
      maxx = a.x[1] > b.x[1] ? a.x[1] : b.x[1],
      miny = a.y[0] < b.y[0] ? a.y[0] : b.y[0],
      maxy = a.y[1] > b.y[1] ? a.y[1] : b.y[1],
      minz = a.z[0] < b.z[0] ? a.z[0] : b.z[0],
      maxz = a.z[1] > b.z[1] ? a.z[1] : b.z[1]
   )
   [ [minx,maxx],[miny,maxy],[minz,maxz] ];

function dm_center(m,center=[true,true,true]) = 
   let(b=dm_bounds(m))
   dm_translate(m,[center.x?(-b.x[0]-b.x[1])/2:0,center.y?(-b.y[0]-b.y[1])/2:0,center.z?(-b.z[0]-b.z[1])/2:0]);
   

// ---------------------------------------------------------------------------------------------------------------------

module dp_polygon(p) {
   polygon(p);
}

module dm_polyhedron(m) {
   polyhedron(m[0],m[1]);
}

// ---------------------------------------------------------------------------------------------------------------------
// == helpers

function lerp(a,b,f) = (1-f)*a + f*b;
function lerp2(a,b,f) = [lerp(a.x,b.x),lerp(a.y,b.y)];
function lerp3(a,b,f) = [lerp(a.x,b.x),lerp(a.y,b.y),lerp(a.z,b.z)];

function _cm_t(t,p0,p1) = 
   let(
      a = abs(p1-p0),
      b = pow(a,0.5),
      c = pow(b,0.5)
   )
   c+t;

function _cm_t2(t,p0,p1) = 
   let(
      a = sqrt(pow(p1.x-p0.x,2)+pow(p1.y-p0.y,2)),
      b = pow(a,0.5),
      c = pow(b,0.5)
   )
   c+t;
   
function _cm_t3(t,p0,p1) = 
   let(
      a = sqrt(pow(p1.x-p0.x,2)+pow(p1.y-p0.y,2)+pow(p1.z-p0-z,2)),
      b = pow(a,0.5),
      c = pow(b,0.5)
   )
   c+t;
   
function dpf_cm_smooth(a,t) = 
   let(
      t0 = 0,
      t1 = _cm_t(t0,a[0],a[1]),
      t2 = _cm_t(t1,a[1],a[2]),
      t3 = _cm_t(t2,a[2],a[3]),
      a1 = (t1-t)/(t1-t0)*a[0] + (t-t0)/(t1-t0)*a[1],
      a2 = (t2-t)/(t2-t1)*a[1] + (t-t1)/(t2-t1)*a[2],
      a3 = (t3-t)/(t3-t2)*a[2] + (t-t2)/(t3-t2)*a[3],
      b1 = (t2-t)/(t2-t0)*a1 + (t-t0)/(t2-t0)*a2,
      b2 = (t3-t)/(t3-t1)*a2 + (t-t1)/(t3-t1)*a3,
      c  = (t2-t)/(t2-t1)*b1 + (t-t1)/(t2-t1)*b2
   )
   c;
   
function dpf_cm_smooth2(a,t) = 
   let(
      t0 = 0,
      t1 = _cm_t2(t0,a[0],a[1]),
      t2 = _cm_t2(t1,a[1],a[2]),
      t3 = _cm_t2(t2,a[2],a[3]),
      a1 = (t1-t)/(t1-t0)*a[0] + (t-t0)/(t1-t0)*a[1],
      a2 = (t2-t)/(t2-t1)*a[1] + (t-t1)/(t2-t1)*a[2],
      a3 = (t3-t)/(t3-t2)*a[2] + (t-t2)/(t3-t2)*a[3],
      b1 = (t2-t)/(t2-t0)*a1 + (t-t0)/(t2-t0)*a2,
      b2 = (t3-t)/(t3-t1)*a2 + (t-t1)/(t3-t1)*a3,
      c  = (t2-t)/(t2-t1)*b1 + (t-t1)/(t2-t1)*b2
   )
   c;
   
function dpf_cm_smooth_alt(a,t) = 
   0.5 * ((2 * a[1]) +
   (-a[0] + a[2]) * t +
   (2*a[0] - 5*a[1] + 4*a[2] - a[3]) * pow(t,2) +
   (-a[0] + 3*a[1]- 3*a[2] + a[3]) * pow(t,3));

function dpf_cm_smooth2(a,f) = [
   dpf_cm_smooth([a[0].x,a[1].x,a[2].x,a[3].x],f),
   dpf_cm_smooth([a[0].y,a[1].y,a[2].y,a[3].y],f)
];

function dpf_cm_smooth3(a,f) = [
   dpf_cm_smooth([a[0].x,a[1].x,a[2].x,a[3].x],f),
   dpf_cm_smooth([a[0].y,a[1].y,a[2].y,a[3].y],f),
   dpf_cm_smooth([a[0].z,a[1].z,a[2].z,a[3].z],f)
];
   
function rand(s) = s!=undef ? rands(0,1,1,s)[0] : rands(0,1,1)[0];

function noise(a,f=1) = rand(a*f);

function noise2(a,f=1) = Coldnoise([a[0]*f,a[1]*f,0]);
function noise3(a,f=1) = Coldnoise([a[0]*f,a[1]*f,a[2]*f]);

// ---------------------------------------------------------------------------------------------------------------------
// third party code snippets

// Catmull-Rom smoothing code by Kit Wallace
cm_weight = [-1, 9, 9, -1] / 16; 

function _cm_interpolate3(arr,n,index_) = [
        arr[(index_ +n - 1) %n].x * cm_weight[0] + 
        arr[index_].x             * cm_weight[1] + 
        arr[(index_ + 1) %n].x    * cm_weight[2] + 
        arr[(index_ + 2) %n].x    * cm_weight[3] ,
        arr[(index_ +n - 1) %n].y * cm_weight[0] + 
        arr[index_].y             * cm_weight[1] + 
        arr[(index_ + 1) %n].y    * cm_weight[2] + 
        arr[(index_ + 2) %n].y    * cm_weight[3] , 
        arr[(index_ +n - 1) %n].z * cm_weight[0] + 
        arr[index_].z             * cm_weight[1] + 
        arr[(index_ + 1) %n].z    * cm_weight[2] + 
        arr[(index_ + 2) %n].z    * cm_weight[3] 
   ]; 

function _cm_subdivide3(arr,index_=0) = 
   index_ >= len(arr) ? [] : 
   concat(
      [arr[index_]], 
      [_cm_interpolate3(arr,len(arr),index_)], 
      _cm_subdivide3(arr,index_+1)
   ); 

function cm_smooth3(arr,depth=4) = 
    depth == 0 ? arr : cm_smooth3(_cm_subdivide3(arr),depth-1); 

// from https://openscadsnippetpad.blogspot.com/2017/05/simple-3d-noise.html
   function Coldnoise(v) =
    let (
      xseed = round(rnd(1e8, -1e8, round(v.x * 1e6))),
      yseed = round(rnd(1e8, -1e8, xseed + round(v.y * 1e6))),
      zseed = round(rnd(1e8, -1e8, yseed + round(v.z * 1e6))),
      noise  =  (rnd(0, 1e8, zseed))%1)
    noise;
   
    function rnd(a = 1, b = 0, s = []) = 
     s == [] ? 
      (rands(min(a, b), max(   a, b), 1)[0]) 
     : 
      (rands(min(a, b), max(a, b), 1, s)[0])
     ; 
  
// from https://en.wikibooks.org/wiki/OpenSCAD_User_Manual/Tips_and_Tricks
   function rotate_from_to(a,b) =      
       let( axis = unit(cross(a,b)) )
       axis*axis >= 0.99 ? 
           transpose([unit(b), axis, cross(axis, unit(b))]) * 
               [unit(a), axis, cross(axis, unit(a))] : 
           identity(3);
   function identity(n) = [for(i=[0:n-1]) [for(j=[0:n-1]) i==j ? 1 : 0] ];
   function unit(v) = norm(v)>0 ? v/norm(v) : undef; 
   function transpose(m) = // m is any rectangular matrix of objects
     [ for(j=[0:len(m[0])-1]) [ for(i=[0:len(m)-1]) m[i][j] ] ];

// from http://kitwallace.tumblr.com/tagged/conway by Kit Wallace, CC BY SA
function vec3(v) = [v.x, v.y, v.z];
function mt_translate(v) = [ [1, 0, 0, 0],
                             [0, 1, 0, 0],
                             [0, 0, 1, 0],
                             [v.x, v.y, v.z, 1  ] ];
function mt_rotate(v) = [ [1,  0,         0,        0],
                          [0,  cos(v.x),  sin(v.x), 0],
                          [0, -sin(v.x),  cos(v.x), 0],
                          [0,  0,         0,        1] ]
                      * [ [ cos(v.y), 0,  -sin(v.y), 0],
                          [0,         1,  0,        0],
                          [ sin(v.y), 0,  cos(v.y), 0],
                          [0,         0,  0,        1] ]
                      * [ [ cos(v.z),  sin(v.z), 0, 0],
                          [-sin(v.z),  cos(v.z), 0, 0],
                          [ 0,         0,        1, 0],
                          [ 0,         0,        0, 1] ];
function mt_transform(v, m)  = vec3([v.x, v.y, v.z, 1] * m);
                            
