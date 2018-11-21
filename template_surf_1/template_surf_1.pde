/*
Managing polyhedral surfaces using half-edge, list based discrete surface representation 
 
 v1
 - defines Face and Surface class
 - imports PLY files 
 - draws a surface by listing all faces and calling drawFace 
 - drawFace draws a face using Processing3D builtin routines (Shape)
 - mousewheel defines camera zoom, mouse motion defines camera rotation
 
 Pascal Romon 2018 
 */


int nVmax = 4000;
int nFmax = 4000;
float distance = 2;  // for camera
Surface S;

void setup() {
  size(500, 500, P3D);
  S = new Surface("cube.ply");  // file in 'data' subfolder
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
  //println(S.getOrderedNeighbours(0));
  S.harmonicFlow();

}

void mouseWheel(MouseEvent event) {  // for zooming in and out
  float e = event.getCount();
  distance += e/100;
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
}

class Surface {
  int nV, nE, nF, nP;   // number of vertices, edges, faces and piramides
  ArrayList<PVector> positions = new ArrayList<PVector>();
  ArrayList<Face> faces = new ArrayList<Face>();
  boolean[][] incidenceVF = new boolean[nVmax][nFmax]; // true iff v is in f 
  boolean[][] adjacency = new boolean[nVmax][nVmax]; // true iff v1 ~ v2 
  ArrayList<IntList> pyramids = new ArrayList<IntList>();

  void drawSurface() {
    int i;
    for (i=0; i<this.nF; i++) {
      drawFace(this, i, 250, 137);
    }
  }

  Surface(String filename) {    // input form external PLY file, has ot be in data folder
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
        } else if(keywords[1].equals("pyramid")) {
          this.nP = int(keywords[2]);
        }
      } else if (keywords[0].equals("end_header")) {
        end_header = true;
      }
      i++;
    }

    // Vertex' 3D coordinates
    for (int j = 0; j < this.nV; j++) {
      String[] keywords = split(lines[i], ' ');
      this.positions.add(new PVector(scalingFactor*float(keywords[0]), 
        scalingFactor*float(keywords[1]), scalingFactor*float(keywords[2])));
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
    
    // pyramid's indexes 
    for (int j=0; j< this.nP; j++) {  // j is the pyramid index
      String[] keywords = split(lines[i], ' ');
      IntList indexes = new IntList();
      int pyramidSize = 4;
      for (int k = 0; k < pyramidSize; k++) {
        int vIndex = int(keywords[k]);
        indexes.append(vIndex);
      }
      pyramids.add(indexes); 
      i++;
    }
    //println(pyramids);
  }
  
  IntList getNeighbours(int vertex) {
    IntList neighbours = new IntList();
    for (int i = 0; i < nV; ++i) {
      if (adjacency[vertex][i]) {
         neighbours.append(i); 
      }
    }
    return neighbours;
  }
  
  boolean getFaceOrientation(int vertex, int neighb1, int neighb2) {
    for (int i = 0; i < faces.size(); ++i) {
      IntList face = faces.get(i).vertices;
      int n = face.size();
      for (int j = 0; j < n; ++j) {
        if (face.get((j - 2 + n)  % n) == vertex && face.get((j - 1 + n) % n) == neighb1  && face.get(j) == neighb2) {
          return true;
        }
      }
    } 
    return false;
  }
  
  IntList getOrderedNeighbours(int vertex) {
    IntList neighbours = getNeighbours(vertex);
    IntList orderedNeighbours = new IntList();
    
    int cur_ind = 0;
    orderedNeighbours.append(neighbours.get(cur_ind));
    neighbours.remove(cur_ind);
    
    while (neighbours.size() > 0) {
      for (int i = 0; i < neighbours.size(); ++i) {
        if (getFaceOrientation(orderedNeighbours.get(cur_ind), vertex, neighbours.get(i))) {
          orderedNeighbours.append(neighbours.get(i));
          neighbours.remove(i);
          ++cur_ind;
          break;
        }
      }
    }
    return orderedNeighbours;
  }
  
  PVector gravityCenter(int vertex) {
    IntList neighbours = getNeighbours(vertex);
    PVector gc = new PVector();
    for(int i = 0; i < neighbours.size(); i++) {
        gc.add(positions.get(neighbours.get(i)));
    }
    if(neighbours.size() == 0) {
      println("No neigbours. Impossible!!!");
      return new PVector();
    }
    return gc.div(neighbours.size());
  }
  
  void harmonicFlow() {
    float tau = 0.001;
    PVector flow = new PVector();
    PVector p = new PVector();
    PVector q = new PVector();
    for(int i = 0; i < nV; i++) {
      p = positions.get(i);
      q = gravityCenter(i);
      flow = PVector.add(p, PVector.mult(PVector.sub(q, p), tau));
      positions.set(i, new PVector(flow.x, flow.y, flow.z));
    }
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
    //println(p);
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
