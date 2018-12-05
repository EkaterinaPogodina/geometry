/*
Managing polyhedral surfaces using list based discrete surface representation 
 
 v1
 - defines Face and Surface class
 - imports PLY files 
 - draws a surface by listing all faces and calling drawFace 
 - drawFace draws a face using Processing3D builtin routines (Shape)
 - mousewheel defines camera zoom, mouse motion defines camera rotation
 
 v2
 - add nextVertex and previousVertex method in Face class
 - vertexDegree method in Surface class
 - orientedNeighbors to find the oriented list of neighbors
 - visualizeNeighbors runs through all pointsvertices and draws the list of neighbors, ordered
 - volume method ; volume is preserved by scaling if scalingVolume == true
 - flow method for arbitrary vector field, with timestep
 - centroid (center of mass)
 - homothety
 - harmonic (combinatorial and probabilistic) vector field
 - to start the flow : 'c' = combinatorial and 'p' = probabilistic
 
 Pascal Romon 2018 
 */


int nVmax = 4000;
int nFmax = 4000;
float distance = 2;  // for camera
Surface S;

int currentVertex = 0;
int currentNeighborIndex = 0;
int timeBetweenVertices = 500; // ms
boolean flow = false;
float tau = 0.05;
char type; // flow type
boolean scalingVolume = true;

void setup() {
  size(500, 500, P3D);
  frameRate(25);
  S = new Surface("mug.ply");  // file in 'data' subfolder
  //S = new Surface("cube.ply");  // file in 'data' subfolder

  /*
  for (int p=0; p < S.nV; p++) {
   //println("p", p, "neighbors", S.incidenceVF[p][0], "next", S.faces.get(0).nextVertex(p));
   IntList neighbors = S.orientedNeighbors(p);
   println("p=", p, neighbors);
   }
   */
}

void draw() {
  background(50);
  camera(distance*width/2.0, distance*height/2.0, distance*(height/2.0) / tan(PI*30.0 / 180.0), 
    width/2.0, height/2.0, 0, 0, 1, 0);
  translate(width/2, height/2, 0);
  rotateX(TWO_PI*mouseY/width);
  rotateY(TWO_PI*mouseX/width);
  // axes
  line(0, 0, 0, width/2, 0, 0);
  line(0, 0, 0, 0, width/2, 0);
  line(0, 0, 0, 0, 0, width/2);

  S.drawSurface();
  if (!flow) {
    S.visualizeNeighbors(timeBetweenVertices);
  } else {
    float volBefore = S.volume();
    println(volBefore*pow(width/2,-1.0/3));  // to compensate scaling from PLY
    switch(type) {
    case 'm': 
      //S.flow(S.MCF(), tau);
      break;
    case 'q':
      S.flow(S.CotanVector(), tau);
      //break;
    case 't':
      S.RMC();
      //break;
    default:
      S.flow(S.harmonic(type), tau);
      break;
    }
    if (scalingVolume) S.homothety(S.centroid(),pow(volBefore/S.volume(),1/3.0));
  }
  //currentVertex = (currentVertex +1) %S.nV;
  //println("p", p, "iVF", S.incidenceVF[p][0], "next", S.faces.get(0).nextVertex(p));
  //IntList neighbors = S.orientedNeighbors(p);
  //println("p=", p, neighbors);
  //for (int i=0; i < neighbors.size(); i++) S.drawVertex(neighbors.get(i), 151, 20);
}

void mouseWheel(MouseEvent event) {  // for zooming in and out
  float e = event.getCount();
  distance += e/100;
}

void keyReleased() {
  if (key == 'c') {  // harmonic flow with combinatorial laplacian
    flow = true;
    type = key;
  }
  if (key == 'p') {  // harmonic flow with probabilistic laplacian
    flow = true;
    type = key;
  }
  if (key == 'm') { // MCF
    flow = true;
    type = key;
  }
  
  if (key == 'q') {
    flow = true;
    type = key;
  }
  if (key == 't') {
    flow = true;
    type = key;
  }
}

///////// CLASSES

class Face {
  // faces are oriented
  IntList vertices = new IntList();

  Face(IntList list) {
    for (int i=0; i<list.size(); i++) {
      this.vertices.append(list.get(i));
    }
  }

  int nextVertex(int p) {  // yields the index of next vertex p in the face of vertex p, -1 if impossible
    int d = vertices.size();   // number of vertices = degree of the face
    int result = -1;
    for (int i=0; i<d && result == -1; i++) {
      if (this.vertices.get(i) == p) result = this.vertices.get((i+1)%d);
    }
    return(result);
  }

  int previousVertex(int p) {  // yields the index of previous vertex p in the face of vertex p, -1 if impossible
    int d = vertices.size();   // number of vertices = degree of the face
    int result = -1;
    for (int i=0; i<d && result == -1; i++) {
      if (this.vertices.get(i) == p) result = this.vertices.get((i-1+d)%d);
    }
    return(result);
  }
}

class Surface {
  int nV, nE, nF;   // number of vertices, edges and faces
  ArrayList<PVector> positions = new ArrayList<PVector>();
  ArrayList<Face> faces = new ArrayList<Face>();
  boolean[][] incidenceVF = new boolean[nVmax][nFmax]; // true iff v is in f 
  boolean[][] adjacency = new boolean[nVmax][nVmax]; // true iff v1 ~ v2 

  Surface(String filename) {    // input form external PLY file, has to be in data folder
    String[] lines = loadStrings(filename);  

    // cleans the matrices first
    for (int i=0; i<nVmax; i++) {
      for (int j=0; j<nFmax; j++) {
        this.incidenceVF[i][j] = false;
      }
      for (int j=0; j<nVmax; j++) {
        this.adjacency[i][j] = false;
      }
    }

    boolean end_header = false;
    String plyTest = "ply";
    if (!lines[0].equals(plyTest)) exit();
    float scalingFactor = width/2;  // assumes coordinates in the PLY file are in [-1,1] 

    int i = 0; // line currently read in the PLY file
    while (!end_header) {
      String[] keywords = split(lines[i], ' ');
      if (keywords[0].equals("element")) {
        if (keywords[1].equals("vertex")) {
          this.nV = int(keywords[2]);
        } else if (keywords[1].equals("face")) {
          this.nF = int(keywords[2]);
        }
      } else if (keywords[0].equals("end_header")) {
        end_header = true;
      }
      i++;
    }
    //println("v=", this.nV, " f=", this.nF);

    // Vertex' 3D coordinates
    for (int j = 0; j < this.nV; j++) {
      String[] keywords = split(lines[i], ' ');
      //println("lines[] " + i + " : " + lines[i]);
      //println("keywords: " + keywords.length);
      this.positions.add(new PVector(scalingFactor*float(keywords[0]), 
        scalingFactor*float(keywords[1]), scalingFactor*float(keywords[2])));
      //println("Vertex " + j + " : " + this.positions.get(j));
      i++;    // increase line number
    }

    // faces' indexes 
    for (int j=0; j< this.nF; j++) {  // j is the face index
      String[] keywords = split(lines[i], ' ');
      IntList indexes = new IntList();
      int degree = int(keywords[0]);
      for (int k=1; k<=degree; k++) {
        int vIndex = int(keywords[k]);
        indexes.append(vIndex);
        this.incidenceVF[vIndex][j] = true;  // vIndex is in face j
      }
      // fills the adjacency matrix
      for (int k=0; k<degree; k++) {
        this.adjacency[indexes.get(k)][indexes.get((k+1) % degree)] = true;
        this.adjacency[indexes.get((k+1) % degree)][indexes.get(k)] = true;
        // to make it symmetric
      }
      Face f = new Face(indexes);
      this.faces.add(f);
      i++;
    }
  }
  void drawSurface() {
    int i;
    for (i=0; i<this.nF; i++) {
      drawFace(this, i, 250, 137);
    }
  }

  void drawVertex(int p, int theStroke, int theFill) {
    if (theStroke == -1) {
      noStroke();
    } else {
      stroke(theStroke);
    }
    if (theFill == -1) {
      //noFill();
      stroke(100, 0, 0);
    } else {
      fill(theFill);
    }
    pushMatrix();
    //println(p,positions.get(p).x, positions.get(p).y, positions.get(p).z);
    translate(positions.get(p).x, positions.get(p).y, positions.get(p).z);
    sphere(10);
    popMatrix();
  }

  int vertexDegree(int p) {      // number of neighbors
    int n = 0;
    for (int q=0; q<nV; q++) {
      if (this.adjacency[p][q]) n++;
    }
    return(n);
  }

  IntList orientedNeighbors(int p) {  // ordered list of neighbors
    // problems when the vertex lies on the boundary
    IntList result = new IntList();

    // get the number of neighbors from adjacency matrix
    int n = this.vertexDegree(p);
    int q = 0; // future first neighbor
    int f; // face containing p and q
    for (q=0; q<nVmax && !adjacency[p][q]; q++);
    result.append(q);
    //println("first neighbor found =", q);

    f = 0;
    while (incidenceVF[p][f] == false || this.faces.get(f).nextVertex(p) != q) f++;
    //println("face found=", f);

    // iterate
    for (int i=1; i<n; i++) {
      // find the next couple vertex-face
      q = this.faces.get(f).previousVertex(p);
      result.append(q);
      for (f=0; f<10 && (incidenceVF[p][f] == false) || (this.faces.get(f).nextVertex(p) != q); f++) {
        //println("B f", f, "iVF", incidenceVF[p][f], "next", faces.get(p).nextVertex(p), "q", q);
      }
    }
    return(result);
  }

  void visualizeNeighbors(int time) {
    // time delay between neighbors (in ms)
    this.drawVertex(currentVertex, 100, 50);
    int q = this.orientedNeighbors(currentVertex).get(currentNeighborIndex);
    this.drawVertex(q, 50, -1);
    delay(time);
    currentNeighborIndex +=1;
    if (currentNeighborIndex == this.vertexDegree(currentVertex)) {
      currentNeighborIndex = 0;
      currentVertex = (currentVertex + 1) %this.nV;
    }
  }

  float volume() {
    float v = 0;
    for (int fIndex=0; fIndex< this.nF; fIndex++) {  // run over all faces
      Face f = faces.get(fIndex);
      PVector p = this.positions.get(f.vertices.get(0)).copy();
      for (int j=1; j < f.vertices.size()-1; j++) {
        PVector q = this.positions.get(f.vertices.get(j));
        PVector r = this.positions.get(f.vertices.get(j+1));
        v += r.dot(p.cross(q));
      }
    }
    return(v/6);
  }

  void flow(ArrayList<PVector> V, float timestep) { // flows the surface along V
    for (int p=0; p < this.nV; p++) 
      this.positions.get(p).add(V.get(p).mult(timestep));
  }
  
  PVector centroid() {    // center of mass of the (vertices of the) surface
    PVector g = new PVector(0,0,0);
    for (int p=0; p < this.nV; p++) g.add(this.positions.get(p));
    return(g.div(this.nV));
  }
  
  void homothety(PVector center, float factor) {
    PVector extra = center.copy().mult(1.0-factor);
    for (int p=0; p < this.nV; p++) 
      this.positions.get(p).mult(factor).add(extra);
  }

  ArrayList<PVector> harmonic(char whichType) { // combinatorial/probabilistic laplacian
    ArrayList<PVector> result = new ArrayList<PVector>();
    for (int p=0; p < this.nV; p++) {
      PVector Vp = new PVector(0, 0, 0);
      for (int q=0; q < this.nV; q++) {
        if (this.adjacency[p][q]) Vp.add(this.positions.get(q)).sub(positions.get(p));
      }
      if (whichType == 'p') Vp.div(this.vertexDegree(p));
      result.add(Vp);
    }
    return(result);
  }

  ArrayList<PVector> CotanVector() {
    ArrayList<PVector> result = new ArrayList<PVector>();
    
    for (int p=0; p < this.nV; p++) {
      PVector Vp = new PVector(0, 0, 0);

      IntList neighbors = orientedNeighbors(p);
      int n = neighbors.size();
      for (int q = 0; q < n; q++) {
        PVector cur = this.positions.get(q);
        PVector next = this.positions.get((q + 1) % n);
        PVector prev = this.positions.get((q - 1 + n) % n);
        PVector point = this.positions.get(q);
    
        PVector pq = PVector.sub(cur, point);
        float psi1 = 1/ tan(arcAngle(pq, PVector.sub(next, cur)));
        float psi2 = 1 / tan(arcAngle(PVector.sub(cur, prev), pq));
        
        Vp.add(PVector.mult(pq, psi2 + psi1));
      }
      result.add(Vp);
    }
 
    return result;
  }

  /*
    Calculates the volume gradient for point p
  */
  ArrayList<PVector> volumeGrad() {
    ArrayList<PVector> result = new ArrayList<PVector>();
    
    for (int p=0; p < this.nV; p++) {

      IntList neighbours = orientedNeighbors(p);
      PVector grad = new PVector();
      PVector pVect = positions.get(p);
      PVector cross = new PVector();
      int n = neighbours.size() - 1; // number of neighbours
      
      for(int i = 0; i < n; i++) {
        PVector.cross(PVector.sub(positions.get(i), pVect), PVector.sub(positions.get((i + 1) % n), pVect), cross);
        grad.add(cross);
      }
      
      result.add(grad.div(6));
    }
 
    return result;
  }
  
  ArrayList<PVector> RMC() {  // renormalized mean curvature
    ArrayList<PVector> H = CotanVector();
    ArrayList<PVector> gV = volumeGrad();
    float lambda = - ndot(H, gV) / ndot(gV, gV);
    return(nsum(H, nmult(gV, lambda)));
  }

}

//////   SUBS

void drawFace(Surface S, int faceIndex, int theStroke, int theFill) {
  // draws face # faceIndex in surface S
  if (theStroke == -1) {
    noStroke();
  } else {
    stroke(theStroke);
  }
  if (theFill == -1) {
    noFill();
  } else {
    fill(theFill);
  }
  Face f = S.faces.get(faceIndex);
  int n = f.vertices.size();
  PVector p;
  beginShape();
  //println(f);
  for (int i=0; i<n; i++) {
    int vertexIndex = f.vertices.get(i);
    p = S.positions.get(vertexIndex);
    vertex(p.x, p.y, p.z);
  }
  endShape();
}

float arcAngle(PVector v1, PVector v2) {
  float r = atan2(v2.y, v2.x) - atan2(v1.y, v1.x);
  if (r<-PI) {
    r += 2*PI;
  }
  if (r>PI) {
    r -= 2*PI;
  }
  return(r);
}

float ndot(ArrayList<PVector> U, ArrayList<PVector> V) {  // scalar product between n-vectors
  float s = 0;
  for (int i=0; i<min(U.size(), V.size()); i++) {
    s += PVector.dot(U.get(i), V.get(i));
  }
  return(s);
}

ArrayList<PVector> nsum(ArrayList<PVector> U, ArrayList<PVector> V) {  // addition between 2 n-vectors
  ArrayList<PVector> W = new ArrayList<PVector>();
  for (int i=0; i<min(U.size(), V.size()); i++) {
    W.add(PVector.add(U.get(i), V.get(i)));
  }
  return(W);
}

ArrayList<PVector> nmult(ArrayList<PVector> U, float lambda) {  // multiplication of n-vector by a float
  ArrayList<PVector> W = new ArrayList<PVector>();
  for (int i=0; i<U.size(); i++) {
    W.add(PVector.mult(U.get(i), lambda));
  }
  return(W);
}
