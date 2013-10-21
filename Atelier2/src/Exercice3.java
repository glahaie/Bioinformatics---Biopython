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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

public class Exercice3 {
    
    public static void main(String[] args) throws FileNotFoundException, IOException, IllegalAlphabetException, IllegalSymbolException {
       
        String fileName = "taxid151340.acc_lst";
        String line;
        List<String> idList = new ArrayList<String>();  // contient la liste des ids des séquences
        int nbID;                                       // nombre total des séquences
        RichSequence hpvRS = null;
        HashMap<String, Float> allSymbolHash = new HashMap<String, Float>();    // contient la fréquence des nucléotides pour toutes les séquences
        
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));

            while ((line = br.readLine()) != null) {
                idList.add(line);
            } 

        }catch (IOException e) {
            System.err.println("Error: " + e);
        }

        nbID = idList.size();
        
        for (String id : idList){
            String seqFile = id + ".gb";
            File f = new File(seqFile);

            // Telechargement des fichiers
            if(! f.exists()) { 
                DownloadGBFile (id);
            }

            // Lecture des fichiers
            hpvRS = readGBFile(seqFile);  
            SymbolList hpvSL = hpvRS.getInternalSymbolList();
            
            int seqLength = hpvSL.length();

            // Calcul des frequences            
            HashMap<String, Integer> symbolHash = new HashMap<String, Integer>();
            
            for (Iterator i = hpvSL.iterator(); i.hasNext();){
                Symbol sym = (Symbol)i.next();
//                System.out.println(sym.getName());
                
                int count = symbolHash.containsKey(sym.getName()) ? symbolHash.get(sym.getName()) : 0;
                symbolHash.put(sym.getName(),count +1);
                
            }
            
            // Affichage
            System.out.println(id);
            
            Iterator<String> keySetIterator = symbolHash.keySet().iterator();
            
            while(keySetIterator.hasNext()){
                String nd = keySetIterator.next();
                float  freq = (float) symbolHash.get(nd) * 100 / seqLength;                
                float count = allSymbolHash.containsKey(nd) ? allSymbolHash.get(nd) : 0;
                
                allSymbolHash.put(nd, count + freq );                        
                System.out.println(nd + " : " + freq);
                
            }
            System.out.println("");        
        }
        
        // Fréquences des nucléotides pour tout le lot 
        System.out.println("ALL : ");
        
        Iterator<String> keySetIterator = allSymbolHash.keySet().iterator();
            
        while(keySetIterator.hasNext()){
            String nd = keySetIterator.next();
            float  freq = (float) allSymbolHash.get(nd)/nbID;
                                          
            System.out.println(nd + " : " + freq);

        }       
    }
    
    public static void DownloadGBFile (String idSeq) throws FileNotFoundException, IOException{
        RichSequence rs = null;
        GenbankRichSequenceDB grsdb = new GenbankRichSequenceDB();
        String GbFile = idSeq + ".gb";
        
        try{
            // Récupération du fichier de GenBank
            rs = grsdb.getRichSequence(idSeq);
//            System.out.println(rs.getDescription());
            
            // Écriture de la séquence dans un fichier GenBank
            OutputStream gbOut = new FileOutputStream(GbFile);
            RichSequence.IOTools.writeGenbank(gbOut, rs, null);
            
        }
        catch(BioException be){
            be.printStackTrace();
            System.exit(-1);
        }
    }
    
    public static RichSequence readGBFile (String gbFile){
        //http://biojava.org/wiki/BioJava:Cookbook:SeqIO:ReadGESBiojavax
        RichSequence rs = null;
        BufferedReader br = null;
        SimpleNamespace ns = null;

        try{
            br = new BufferedReader(new FileReader(gbFile));
            ns = new SimpleNamespace("biojava");
     
            RichSequenceIterator rsi = RichSequence.IOTools.readGenbankDNA(br,ns);
            rs = rsi.nextRichSequence();
//            System.out.println(rs.getName());
        
        }
        catch(Exception be){
                be.printStackTrace();
                System.exit(-1);
        }      
        return rs ;
    }   
}
