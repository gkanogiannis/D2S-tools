package agkanogiannis.hcluster;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.StringWriter;
import java.io.Writer;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeModel;

public class PhylipWriter {
	protected Writer writer;
    protected Object output;
    
    public PhylipWriter() {}

    public void setOutput(Writer out) {
        writer = out;
        output = null;
    }
    
    public void setOutput(File file) throws IOException {
        setOutput(file, false);
    }
    
    public void setOutput(File file, boolean append) throws IOException {
        if(file == null)
            writer = new BufferedWriter(new OutputStreamWriter(System.out));
        else
            writer = new BufferedWriter(new FileWriter(file.getCanonicalPath(),
                append));
        output = file;
    }

    public void setOutput(OutputStream out) {
        writer = new BufferedWriter(new OutputStreamWriter(out));
        output = null;
    }

    public String toString() {
        return(writer == null ? null : writer.toString());
    }

    public void flush() throws IOException {
        if(writer != null)
            writer.flush();
    }

    public void close() throws IOException {
        if(writer != null)
            writer.close();
    }

    public void open() throws IOException {
        if(output == null)
            return;
        else if(output instanceof File)
            setOutput((File)output);
        else if(writer instanceof StringWriter)
            setOutput(new StringWriter());
    }
    
    public Object write(Object obj) throws IOException {
        try {
            DefaultMutableTreeNode root;
            if(obj instanceof TreeModel)
                root = (DefaultMutableTreeNode)((TreeModel)obj).getRoot();
            else if(obj instanceof TreeDataModel)
                root = ((TreeDataModel)obj).getRoot();
            else
                throw new ClassCastException("Tree must be TreeModel or TreeDataModel: "
                    + (obj == null ? "null" : String.valueOf(obj.getClass())));
            writeTree(root);
            writer.write(";\n");
        }
        finally {
            writer.flush();
            writer.close();
        }
        if(writer instanceof StringWriter)
            return writer.toString();
        else
            return null;
    }

    /**
     * Recursively write tree nodes.
     */
    protected void writeTree(DefaultMutableTreeNode node) throws IOException {
        if(!node.isLeaf())
            writer.write("(\n");

        //Children
        for(int i = 0, numChildren = node.getChildCount(); i < numChildren; i++) {
            writeTree((Clade)node.getChildAt(i));
            if(i < numChildren - 1)
                writer.write(",\n");
        }
        if(!node.isLeaf())
            writer.write(")\n");

        //Label
        if(!node.isRoot() && node.getUserObject() != null)
            writer.write(node.toString());

        //Branch length
        if(!node.isRoot() && node instanceof Clade) {
            writer.write(":");
            Clade branch = (Clade)node;
            double length = branch.getBranchLength();
            writer.write(String.valueOf(length));
        }
    }
}