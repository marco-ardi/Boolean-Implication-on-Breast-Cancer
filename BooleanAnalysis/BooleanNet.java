import java.util.HashMap;
import java.util.Vector;

public class BooleanNet
{
    public static void main(String[] args) throws Exception
    {
        String inputFile=null;
        String outputFile="implications.txt";
        String gene1=null;
        String gene2=null;
        String type="any";
        double statThresh=6.0;
        double pvalThresh=0.01;

        //Read input parameters
        int i=0;
        for (i=0;i<args.length;i++)
        {
            switch (args[i])
            {
                case "-i" -> inputFile = args[++i];
                case "-g1" -> gene1 = args[++i];
                case "-g2" -> gene2 = args[++i];
                case "-rt" -> type = args[++i];
                case "-s" -> statThresh = Double.parseDouble(args[++i]);
                case "-p" -> pvalThresh = Double.parseDouble(args[++i]);
                case "-o" -> outputFile = args[++i];
                default -> {
                    System.out.println("Error! Unrecognizable command '" + args[i]);
                    printHelp();
                    System.exit(1);
                }
            }
        }

        //Error in case network file is missing
        if(inputFile==null)
        {
            System.out.println("Error! No expression matrix has been specified!\n");
            printHelp();
            System.exit(1);
        }

        //Read expression data
        FileManager fm=new FileManager();
        HashMap<String,int[]> mapExpressions=fm.readDiscretizedExpressionFile(inputFile);

        //Start boolean analysis
        BooleanAnalyzer sa=new BooleanAnalyzer();
        fm.initImplicationFile(outputFile);

        //Find implications and write results
        int numTotalImplications=0;
        if(gene1==null && gene2==null)
        {
            int numProcessedGenes=0;
            int numTotalGenes= mapExpressions.size();
            int numChecks=10;
            int[] checkPoints=new int[numChecks+1];
            for(i=0;i<checkPoints.length-1;i++)
                checkPoints[i]=numTotalGenes/(checkPoints.length-1)*i;
            checkPoints[i]=numTotalGenes;
            int currentIndexCheck=0;
            for(String gene : mapExpressions.keySet())
            {
                if(numProcessedGenes==checkPoints[currentIndexCheck])
                {
                    System.out.println("Finding implications of "+type+" type for all genes..."+currentIndexCheck*numChecks+"%");
                    currentIndexCheck++;
                }
                Vector<Vector<String>> listImplications=sa.getImplications(gene,type,mapExpressions,statThresh,pvalThresh);
                fm.writeImplications(listImplications,outputFile);
                numTotalImplications+=listImplications.size();
                numProcessedGenes++;
            }
        }
        else if(gene2==null)
        {
            if(mapExpressions.containsKey(gene1))
            {
                System.out.println("Finding implications of "+type+" type for gene "+gene1+"...");
                Vector<Vector<String>> listImplications=sa.getImplications(gene1,type,mapExpressions,statThresh,pvalThresh);
                fm.writeImplications(listImplications,outputFile);
                numTotalImplications+=listImplications.size();
            }
            else
                System.out.println("Warning! Gene "+gene1+" is not present in input data. No implications for "+gene1+" have been retrieved");
        }
        else
        {
            if(mapExpressions.containsKey(gene1) && mapExpressions.containsKey(gene2))
            {
                System.out.println("Finding implications of "+type+" type between gene "+gene1+" and gene "+gene2+"...");
                Vector<Vector<String>> listImplications=sa.getImplications(gene1,gene2,type,mapExpressions,statThresh,pvalThresh);
                fm.writeImplications(listImplications,outputFile);
                numTotalImplications+=listImplications.size();
            }
            else
            {
                if(!mapExpressions.containsKey(gene1))
                    System.out.println("Warning! Gene "+gene1+" is not present in input data. No implications for "+gene1+" have been retrieved");
                if(!mapExpressions.containsKey(gene2))
                    System.out.println("Warning! Gene "+gene2+" is not present in input data. No implications for "+gene2+" have been retrieved");
            }
        }
        System.out.println("Done! Found "+numTotalImplications+" implications.");
    }

    private static void printHelp()
    {
        String help = "Usage: java BooleanNet -i <discretizedMatrixFile> [-g1 <geneA> -g2 <geneB> -rt <relType> ";
        help+="-s <statisticThreshold> -p <pvalThreshold> -o <outputFile>]\n\n";
        help+="REQUIRED PARAMETERS:\n";
        help+="-i\tDiscretized expression matrix file\n\n";
        help+="OPTIONAL PARAMETERS:\n";
        help+="-g1\tFirst gene (default=all genes)\n";
        help+="-g2\tSecond gene (default=all genes)\n";
        help+="-rt\tRelationship type (default='any', possible values='any','low-low','low-high','high-low','high-high',";
        help+="'equivalent','opposite')\n";
        help+="-s\tMinimum value of the statistic S of an implication to be considered significant (default=6.0)\n";
        help+="-p\tMaximum p-value P of an implication to be considered significant (default=0.01)\n";
        help+="-o\tOutput file where implications found will be saved (default=implications.txt)\n";
        System.out.println(help);
    }
}
