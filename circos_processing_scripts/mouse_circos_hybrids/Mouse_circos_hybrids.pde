import processing.pdf.*;
Table genome_info, sp, lz, dip, rs, spSig, lzSig, dipSig, rsSig;
PFont theBest, theBold;
ChromRing def;
GeneRing spRing, lzRing, dipRing, rsRing;
SigRing spSigRing, lzSigRing, dipSigRing, rsSigRing;
int ind = 0;

//int sp_color = #E89C87;
int sp_color = #D85531;
int sp_color_dark = #D85531;
//int lz_color = #F4D5A4;
int lz_color = #ECB45B;
int lz_color_dark = #ECB45B;
//int dip_color = #9AC1AE;
int dip_color = #5D987B;
int dip_color_dark = #5D987B;
//int rs_color = #6C719D;
int rs_color = #3D405B;
int rs_color_dark = #3D405B;

int tick_length = 15;
int tick_spacer = 12;
int tick_centerRing = 225;
float sp_position    =  tick_centerRing + (-1.5*   (tick_spacer + tick_length));
float lz_position    =  tick_centerRing + (0.5*    (tick_spacer + tick_length));
float dip_position   =  tick_centerRing + (2.5*    (tick_spacer + tick_length));
float rs_position    =  tick_centerRing + (4.5*    (tick_spacer + tick_length));
float spSig_position =  tick_centerRing + (-2.3*   (tick_spacer + tick_length));
float lzSig_position =  tick_centerRing + (-0.4*   (tick_spacer + tick_length));
float dipSig_position = tick_centerRing + (1.6*    (tick_spacer + tick_length));
float rsSig_position  = tick_centerRing + (3.6*    (tick_spacer + tick_length));

void setup() {
  size(800, 800);
  strokeCap(SQUARE);
  ellipseMode(CENTER);

  genome_info = loadTable("mmus_chromosomes.csv", "header");

  //beginRecord(PDF, "DE_genes_hybrids_both_parents_logFC_MusP.pdf");
  //sp = loadTable("hybrid_SP.csv", "header");
  //lz = loadTable("hybrid_LZ.csv", "header");
  //dip = loadTable("hybrid_DIP.csv", "header");
  //rs = loadTable("hybrid_RS.csv", "header");
  //spSig = loadTable("Hyb.Comb.SP_sig.csv", "header");
  //lzSig = loadTable("Hyb.Comb.LZ_sig.csv", "header");
  //dipSig = loadTable("Hyb.Comb.DIP_sig.csv", "header");
  //rsSig = loadTable("Hyb.Comb.RS_sig.csv", "header");
  
    beginRecord(PDF, "mouse_DE_genes_hybrids_trans_logFC_MusP.pdf");

  sp = loadTable("hybrid_mouse_trans_SP_DE_Genes.csv", "header");
  lz = loadTable("hybrid_mouse_trans_LZ_DE_Genes.csv", "header");
  dip = loadTable("hybrid_mouse_trans_DIP_DE_Genes.csv", "header");
  rs = loadTable("hybrid_mouse_trans_RS_DE_Genes.csv", "header");
  spSig = loadTable("Mouse.SP_sig.csv", "header");
  lzSig = loadTable("Mouse.LZ_sig.csv", "header");
  dipSig = loadTable("Mouse.DIP_sig.csv", "header");
  rsSig = loadTable("Mouse.RS_sig.csv", "header");

  //noSmooth();
  theBest = createFont("HelveticaNeue-Thin", 16);
  theBold = createFont("Avenir-Heavy", 12);
  textFont(theBest);

  genome_info.addColumn("length", Table.FLOAT);
  float chromNum = genome_info.getRowCount();
  //float genomeLength = 2.098794908E9;
  float genomeLength = 2.725521370E9;
  genomeLength = 2.725521370E9 + genomeLength/100;
  float chromSpacer = genomeLength/200;
  genomeLength = genomeLength + (chromSpacer*chromNum);
  def = new ChromRing(genomeLength, chromSpacer);
  spRing = new GeneRing(genomeLength, chromSpacer, sp_position, sp, tick_length, sp_color_dark, sp_color);
  lzRing = new GeneRing(genomeLength, chromSpacer, lz_position, lz, tick_length, lz_color_dark, lz_color);
  dipRing = new GeneRing(genomeLength, chromSpacer, dip_position, dip, tick_length, dip_color_dark, dip_color);
  rsRing = new GeneRing(genomeLength, chromSpacer, rs_position, rs, tick_length, rs_color_dark, rs_color);
    
  spSigRing = new SigRing(genomeLength, chromSpacer, spSig_position, spSig, sp_color_dark, sp_color);
  lzSigRing = new SigRing(genomeLength, chromSpacer, lzSig_position, lzSig, lz_color_dark, lz_color);
  dipSigRing = new SigRing(genomeLength, chromSpacer, dipSig_position, dipSig, dip_color_dark, dip_color);
  rsSigRing = new SigRing(genomeLength, chromSpacer, rsSig_position, rsSig, rs_color_dark, rs_color);
}

void draw() {
  background(255);
  strokeCap(SQUARE);
  ellipseMode(CENTER);

  spRing.display();
  lzRing.display();
  dipRing.display();
  rsRing.display();

    
  textFont(theBold);
  spSigRing.display();
  lzSigRing.display();
  dipSigRing.display();
  rsSigRing.display();
  textFont(theBest);

  def.display();

  //legend setup
  int legendH = 210;
  int legendW = 220;
  int squareH = 17;
  float spacerY = 80/7;
  float spacerX = 5;
  strokeWeight(1);

  noFill();
  rect(400-(legendW/2), 400-(legendH/2), legendW, legendH);
  ellipseMode(CENTER);
  textAlign(LEFT);
  noStroke();

  //heading arrows and text;
  fill(0);
  //only arrow
  float centerX = 400-(legendW/2) + spacerY + 0.75*squareH/2 + squareH/12;
  float centerY = 400-(legendH/2) + spacerY + 0*(spacerY + squareH);
  //text
  fill(0);
  
  textFont(theBold);
  textSize(16);
  text("DE genes by sample type", 400 - textWidth("DE genes by sample type")/2, 400-(legendH/2) + spacerY + 0*(spacerY + squareH) + textAscent() - 2);
  textFont(theBest);
  //underline
  rect(400-(legendW/2) + 4, 400 - (legendH/2) + spacerY + (spacerY + squareH) - 6, legendW - 8, squareH/12);

  //SP boxes and text; sp_color_dark, sp_color
  fill(sp_color);
  rect(400-(legendW/2) + 1.5*spacerX, 400-(legendH/2) + spacerY + (spacerY + squareH), squareH, squareH);
  fill(0);
  text("Spermatogonia (SP)", 400-(legendW/2) + spacerY + 1*(spacerX + squareH), 400-(legendH/2) + spacerY + (spacerY + squareH) + textAscent() - 2);

  //LZ boxes and text; lz_color_dark, lz_color
  fill(lz_color);
  rect(400-(legendW/2) + 1.5*spacerX, 400-(legendH/2) + spacerY + 2*(spacerY + squareH), squareH, squareH);
  fill(0);
  text("Leptotene/zygotene (LZ)", 400-(legendW/2) + spacerY + 1*(spacerX + squareH), 400-(legendH/2) + spacerY + 2*(spacerY + squareH) + textAscent()-2);

  //DIP boxes and text; #dip_color_dark, dip_color
  fill(dip_color);
  rect(400-(legendW/2) + 1.5*spacerX, 400-(legendH/2) + spacerY + 3*(spacerY + squareH), squareH, squareH);
  fill(0);
  text("Diplotene (DIP)", 400-(legendW/2) + spacerY + 1*(spacerX + squareH), 400-(legendH/2) + spacerY +  3*(spacerY + squareH) + textAscent() - 2);

  //RS boxes and text; #rs_color_dark, rs_color
  fill(rs_color);
  rect(400-(legendW/2) + 1.5*spacerX, 400-(legendH/2) + spacerY + 4*(spacerY + squareH), squareH, squareH);
  fill(0);
  text("Round spermatids (RS)", 400-(legendW/2) + spacerY + 1*(spacerX + squareH), 400-(legendH/2) + spacerY +  4*(spacerY + squareH) + textAscent() - 2);

  //over-enriched and under-enriched
  fill(#000000);
  textFont(theBold);
  textSize(16);
  textAlign(LEFT);
  text("+", 400-(legendW/2) + 1.5*spacerX + 0.5*squareH - 0.5*textWidth("+"), 400-(legendH/2) + spacerY +  5*(spacerY + squareH) + textAscent() - 2);
  text("-", 400-(legendW/2) + 1.5*spacerX + 0.5*squareH - 0.35*textWidth("-"), 400-(legendH/2) + spacerY +  6*(spacerY + squareH) + textAscent() - 2);

  textFont(theBest);
  textSize(16);
  textAlign(LEFT);
  text("Over-enriched for DEs", 400-(legendW/2) + spacerY + 1*(spacerX + squareH), 400-(legendH/2) + spacerY +  5*(spacerY + squareH) + textAscent() - 2);
  text("Under-enriched for DEs", 400-(legendW/2) + spacerY + 1*(spacerX + squareH), 400-(legendH/2) + spacerY +  6*(spacerY + squareH) + textAscent() - 2);


  //ring labelling
  textAlign(CENTER);

  textFont(theBold);
  textSize(14);
  fill(sp_color);
  text("SP", 400 + sp_position, 400);
  fill(lz_color);
  text("LZ", 400 + lz_position, 400);
  fill(dip_color);
  text("DIP", 400 + dip_position, 400);
  fill(rs_color);
  text("RS", 400 + rs_position, 400);

  fill(0);
  text("Chr", 400 + 375, 400);



  endRecord();
  noLoop();
}

class ChromRing {

  float genomeLength;
  float chromSpacer;

  ChromRing(float gen_len, float chrom_space) {
    genomeLength = gen_len;
    chromSpacer = chrom_space;
  }

  void display() {
    float prevEnd = 0;
    String search;
    for (TableRow row : genome_info.rows()) {
      strokeWeight(5);
      //for (TableRow row : genome_info.matchRows("^16$", "chromosome")) {
      search = row.getString("chromosome");
      float start = row.getFloat("start");
      float end = row.getFloat("end");
      row.setFloat("length", end - start);
      float thisLength = row.getFloat("length");

      float chromStart = map(start + prevEnd, 1, genomeLength, 0.03, 2*PI);
      float chromStop = map(end + prevEnd, 1, genomeLength, 0.03, 2*PI);

      if ((ind % 2) == 0) {
        stroke(#E3E3E3);
      } else {
        stroke(#A7A7A7);
      }
      noFill();
      arc(400, 400, 750, 750, chromStart, chromStop, OPEN);
      ind++;

      pushMatrix();
      translate(400, 400);
      stroke(0);
      rotate(-PI/2);
      rotate((chromStart + chromStop) / 2);
      textAlign(CENTER);
      fill(0);
      textSize(12);
      if (search.equals("10") || search.equals("11") || search.equals("12") || search.equals("13") || search.equals("14") || search.equals("15") || search.equals("16") || search.equals("17") || search.equals("18") || search.equals("19") || search.equals("X") || search.equals("Y")) {
        pushMatrix();
        scale(-1, -1);
        text(search, 0, -(365 + 22 + 7) + (textAscent()/2));
        popMatrix();
      } else {
        text(search, 0, 365 + 22 + (textAscent()/2));
      }
      popMatrix();


      if (search.equals("1")) {
        prevEnd = thisLength + chromSpacer;
      } else {
        prevEnd += thisLength + chromSpacer;
      }
    }
  }
}


class GeneRing {

  float genomeLength;
  float chromSpacer;
  float radius;
  Table file;
  int thickneDomP;
  int ringColorUp;
  int ringColorDown;

  GeneRing(float gen_len, float chrom_space, float rad, Table de_file, int thick, int hexUp, int hexDown) {
    genomeLength = gen_len;
    chromSpacer = chrom_space;
    radius = rad;
    file = de_file;
    thickneDomP = thick;
    ringColorUp = hexUp;
    ringColorDown = hexDown;
  }

  void display() {
    noFill();
    float prevEnd = 0;
    String search;
    int currentChr = 0;

    for (TableRow row : genome_info.rows()) {
      //for (TableRow row : genome_info.matchRows("^16$", "chromosome")) {
      search = row.getString("chromosome");
      float start = row.getFloat("start");
      float end = row.getFloat("end");
      row.setFloat("length", end - start);
      float thisLength = row.getFloat("length");
      println(thisLength);

      for (TableRow row2 : file.matchRows("^" + search + "$", "chr")) {

        //for (TableRow row2 : de_genes.matchRows("^2$", "chr")) {
        float localStart = row2.getFloat("start") + prevEnd;
        float localEnd = row2.getFloat("stop") + prevEnd;
        float mappedStart = map(localStart, 1, genomeLength, 0.03, 2*PI);
        float mappedStop = map(localEnd, 1, genomeLength, 0.03, 2*PI);
        println(row2.getString("chr") + " " + row2.getString("start") + " " + mappedStart);

        //up or down
        float fc = row2.getFloat("logFC.x");
        if (fc > 0) {
          stroke(ringColorUp);
        } else if (fc < 0) {
          stroke(ringColorDown);
        }


        //arc(400, 400, (300+radius+thickneDomP), (300+radius+thickneDomP), mappedStart, mappedStop*1.001, OPEN);
        //println(mappedStart + " " + mappedStop);
        //        println((2 + ((localEnd - localStart)/300000))*1.5);
        pushMatrix();
        translate(400, 400);
        rotate(-PI/2);
        rotate(mappedStart);
        //strokeWeight((2 + ((localEnd - localStart)/1500000))*1.5);
        strokeWeight(0.5);
        line(0, radius+thickneDomP, 0, radius-thickneDomP);
        popMatrix();
      }
      if (currentChr == 0) {
        prevEnd = thisLength + chromSpacer;
      } else {
        prevEnd += thisLength + chromSpacer;
      }
      currentChr++;
    }
  }
}

class SigRing {

  float genomeLength;
  float chromSpacer;
  float radius;
  Table sigFile;
  int ringColorUp;
  int ringColorDown;

  SigRing(float gen_len, float chrom_space, float rad, Table sig_file, int hexUp, int hexDown) {
    genomeLength = gen_len;
    chromSpacer = chrom_space;
    radius = rad;
    sigFile = sig_file;
    ringColorUp = hexUp;
    ringColorDown = hexDown;
  }


  void display() {
    float prevEnd = chromSpacer;
    float currentMidLength = 0;
    String search;
    int currentChr = 0;

    for (TableRow row : genome_info.rows()) {
      search = row.getString("chromosome");
      float start = row.getFloat("start");
      float end = row.getFloat("end");
      row.setFloat("length", end - start);
      float thisLength = row.getFloat("length");
      println(thisLength);
      float midLength = thisLength/2;
      currentMidLength = prevEnd + midLength;

      for (TableRow row2 : sigFile.matchRows("^" + search + "$", "chr")) {

        float mappedMid = map(currentMidLength, 1, genomeLength, 0.03, 2*PI);

        //up or down
        String direction = row2.getString("direction");
        if (direction.equals("Over")) {
          fill(ringColorUp);
        } else {
          fill(ringColorDown);
        }

        String chrText = row2.getString("ringText");

        pushMatrix();
        translate(400, 400);
        rotate(-PI/2);
        rotate(mappedMid);
        textAlign(CENTER);
        textSize(12);
        text(chrText, 0, radius);
        popMatrix();
        textSize(16);
      }
      if (currentChr == 0) {
        prevEnd = thisLength + chromSpacer;
      } else {
        prevEnd += thisLength + chromSpacer;
      }
      currentChr++;
    }
  }
}
