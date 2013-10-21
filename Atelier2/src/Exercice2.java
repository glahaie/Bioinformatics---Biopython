/*
 * INF4500 Bioinformatique
 * Automne 2013
 * UQAM
 */

/**
 * @author Mohamed AMine Remita <remita.mohamed_amine@uqam.ca>
 * @author Abdoulaye Baniré Diallo <diallo.abdoulaye@uqam.ca>
 */

import java.io.*;
import java.util.Iterator;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.AlignmentPair;
import org.biojava.bio.alignment.NeedlemanWunsch;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.Note;
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;
import org.biojavax.bio.db.ncbi.GenpeptRichSequenceDB;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.ontology.ComparableTerm;

public class Exercice2 {
    public static void main(String[] args) throws FileNotFoundException, IOException, IllegalAlphabetException {
        
         // 1. Telechargement du fichier NM_001123784.gb
        
        RichSequence rs = null;
        String SeqID = "FR872717";
        String GbFile = SeqID + ".gb";
        
        GenbankRichSequenceDB grsdb = new GenbankRichSequenceDB();
        
        try{
            // Récupération du fichier de GenBank
            rs = grsdb.getRichSequence(SeqID);
            
            // Écriture de la séquence dans un fichier GenBank
            OutputStream gbOut = new FileOutputStream(GbFile);
            RichSequence.IOTools.writeGenbank(gbOut, rs, null);
        }
        catch(BioException be){
            be.printStackTrace();
            System.exit(-1);
        }
        
        // 2. Gènes (gene or CDS)
        // http://biojava.org/wiki/BioJava:Cookbook:Annotations:List2
        
        //Filter the sequence on CDS features
        FeatureFilter ff = new FeatureFilter.ByType("CDS");
        FeatureHolder fh = rs.filter(ff);
        
        System.out.println("#Gene\tlocation\tprotein id\tproduct");
        
        for (Iterator <Feature> i = fh.features(); i.hasNext();){
            RichFeature rf = (RichFeature)i.next();

            //Get the location of the feature
            String featureLocation = rf.getLocation().toString();

            //Get the annotation of the feature
            RichAnnotation ra = (RichAnnotation)rf.getAnnotation();

            //Use BioJava defined ComparableTerms 
            ComparableTerm geneTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("gene");
            ComparableTerm productTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("product");
            ComparableTerm proteinIDTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("protein_id");

            //Create empty strings
            String gene = "";   
            String product = "";   
            String proteinID = "";
            
            //Iterate through the notes in the annotation 
            for (Iterator <Note> it = ra.getNoteSet().iterator(); it.hasNext();){
                Note note = it.next();

            //Check each note to see if it matches one of the required ComparableTerms
                
                if(note.getTerm().equals(productTerm)){
                product = note.getValue().toString();
                }
                if(note.getTerm().equals(geneTerm)){
                gene = note.getValue().toString();
                }
                
                if(note.getTerm().equals(proteinIDTerm)){
                proteinID = note.getValue().toString();
                }
            }
            //Outout the feature information
            System.out.println(gene  +  "\t" + featureLocation + "\t" + proteinID + "\t" + product);
        }
        
        // 3. Extraction et traduction de la séquence du gène L1
        // Le gène L1 est situé dans cette position 5771..7276
        
        SymbolList cds = rs.getInternalSymbolList().subList(5771, 7276);
        // transcription en ARN
        cds = DNATools.toRNA(cds);
        // traduction
        SymbolList protSL = RNATools.translate(cds);
        
//        System.out.println("\n" + protSL.seqString());
       
        // 4. Telechargement de séquence de la proteine L1
        // Le numéro d'accession de la protine L1 est CCB84764
        RichSequence protRS = null;
        try{
            // Récupération du fichier de la base NCBI PROTEIN
            GenpeptRichSequenceDB prsdb = new GenpeptRichSequenceDB();          // different de l'objet GenbankRichSequenceDB
            
            protRS = prsdb.getRichSequence("CCB84764.1");
            
            GbFile = "CCB84764.fa";
            // Écriture de la séquence dans un fichier Fasta
            OutputStream gbOut = new FileOutputStream(GbFile);
            RichSequence.IOTools.writeFasta(gbOut, protRS, null);
        }
        catch(BioException be){
            be.printStackTrace();
            System.exit(-1);
        }
      
        // 5. Comparaison => Alignement global de 2 séquences 
        // http://biojava.org/wiki/BioJava:CookBook:DP:PairWise2
        
        try {
            // The alphabet of the sequences. For this example DNA is choosen.
            FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
            
            // Read the substitution matrix file. 
            String matrixFile = "BLOSUM40.mtx";
            SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(matrixFile));
   
            // Define the default costs for sequence manipulation for the global alignment.
//            SequenceAlignment aligner = new NeedlemanWunsch(          // La classe SequenceAlignment est absente de Biojava 1.8
            NeedlemanWunsch aligner = new NeedlemanWunsch( 
                (short) 0, 	// match
                (short) 3,	// replace
                (short) 2,      // insert
                (short) 2,	// delete
                (short) 1,      // gapExtend
                matrix 	// SubstitutionMatrix
            );
            
            // Perform an alignment and save the results.
            
            SymbolList protSL2 = protRS.getInternalSymbolList();
                              
            protSL2 = protSL2.subList(1, protSL2.seqString().length()-1);
              
            AlignmentPair alignPair = aligner.pairwiseAlignment(
                protSL, // first sequence
                protSL2 // second one
            );
            
             //Print the alignment to the screen
            System.out.println("\nGlobal alignment with Needleman-Wunsch:\n" + alignPair.formatOutput());
        
//            System.out.println(protSL.getAlphabet().getName());
//            System.out.println(protSL2.getAlphabet().getName());
//            System.out.println(matrix.getAlphabet().getName());
            
        } catch (Exception exc) {
            exc.printStackTrace();
        }
                
    }
}
