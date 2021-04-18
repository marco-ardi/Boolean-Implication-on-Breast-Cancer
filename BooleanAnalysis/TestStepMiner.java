import java.util.Vector;

public class TestStepMiner
{
    public static void main(String[] args) throws Exception
    {
        //String expressionFile="Script/Example/tcga.txt";
        String expressionFile="Data/bladder_expr.txt";
        //String discretizedMatrixFile="Output/tcga_disc.txt";
        String discretizedMatrixFile="Output/bladder_expr_disc.txt";
        double gap=0.5;

        FileManager fm=new FileManager();

        //Read expression data
        Vector<double[]> mapExpressions=fm.readExpressionFile(expressionFile);

        //Start boolean analysis
        BooleanAnalyzer sa=new BooleanAnalyzer();

        //Run StepMiner to discretize expression values
        Vector<int[]> listDiscretizedValues=sa.discretizeMatrix(mapExpressions,gap);

        //Write output results to file
        fm.writeDiscretizedMatrix(listDiscretizedValues,discretizedMatrixFile);

    }
}
