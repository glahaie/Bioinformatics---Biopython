/*
 * INF4500 Bioinformatique
 * Automne 2013
 * UQAM
 */

/**
 * @author Mohamed AMine Remita <remita.mohamed_amine@uqam.ca>
 * @author Abdoulaye Baniré Diallo <diallo.abdoulaye@uqam.ca>
 */

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;
import org.biojavax.bio.seq.RichSequence;


public class Exercice1 {
   
    public static void main(String[] args) throws IOException, IllegalAlphabetException {
        
        // 1. Telechargement du fichier NM_001123784.gb
        // http://biojava.org/wiki/BioJava:CookBook:ExternalSources:NCBIFetch
        
        RichSequence rs = null;
        String SeqID = "NM_001123784";
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
         
        
         // 2. Affichage des informations 
        
        String SeqInfo = "Numéro d'accession  : " + rs.getName()+ " \n" +
                         "GI : " +  rs.getIdentifier() + " \n" +
                         "Version : " +  rs.getVersion() + " \n"+
                         "Nom : " + rs.getDescription() + " \n" +
                         "Taille : " + rs.length() + " \n" +
                         "Taxon : " + rs.getTaxon();
        
        System.out.println(SeqInfo);
        
        
        // 3. Traduction de la séquence
        // http://biojava.org/wiki/BioJava:Cookbook:Translation
        // Traduction de toute la séquence
        
        SymbolList sl = rs.getInternalSymbolList();
        // transcription en ARN
        sl = DNATools.toRNA(sl);
        // traduction
        sl = RNATools.translate(sl);
        
        System.out.println("\n" + sl.seqString());
        // On trouve des (*) dans cette séquence qui correspondent à des codons Stop
        // ce qui indique que ce n'est pas la totalité de cet ARNm sera traduite,
        // Pour cela on essaye de traduire seulement le CDS.
        
        // Traduction du CDS seulement qui se situe entre 452 et 3496
          
        SymbolList cds = rs.getInternalSymbolList().subList(452, 3496);
        // transcription en ARN
        cds = DNATools.toRNA(cds);
        // traduction
        SymbolList protSL = RNATools.translate(cds);
        System.out.println("\n" + protSL.seqString());
        
        
        // 4. Enregistrement de la séquence de la PROTEINE
        // http://biojava.org/wiki/BioJava:Cookbook:SeqIO:WriteInFasta
        
        String ProtFile = SeqID + "_prot.fa";
        OutputStream faOut = new FileOutputStream(ProtFile);
        
        RichSequence protRS = RichSequence.Tools.createRichSequence(rs.getName() + " protein", protSL);
        RichSequence.IOTools.writeFasta(faOut, protRS, null);
    }
}
