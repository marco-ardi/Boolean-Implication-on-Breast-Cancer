import java.util.HashMap;
import java.util.Vector;

public class TestBooleanNet
{
    public static void main(String[] args) throws Exception
    {
        String expressionFile="Output/bladder_expr_disc.txt";
        String implicationFile="Output/bladder_expr_impl.txt";
        String gene1=null;
        String gene2=null;
        String type="opposite";
        double statThresh=6.0;
        double pvalThresh=0.01;

        FileManager fm=new FileManager();

        //Read expression data
        HashMap<String,int[]> mapExpressions=fm.readDiscretizedExpressionFile(expressionFile);

        //Start boolean analysis
        BooleanAnalyzer sa=new BooleanAnalyzer();
        fm.initImplicationFile(implicationFile);

        //Find implications and write results
        int numTotalImplications=0;
        if(gene1==null && gene2==null)
        {
            int numProcessedGenes=0;
            int numTotalGenes= mapExpressions.size();
            int numChecks=10;
            int[] checkPoints=new int[numChecks+1];
            int i=0;
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
                fm.writeImplications(listImplications,implicationFile);
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
                fm.writeImplications(listImplications,implicationFile);
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
                fm.writeImplications(listImplications,implicationFile);
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
}
