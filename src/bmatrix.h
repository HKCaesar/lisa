#ifndef BMATRIX_H
#define BMATRIX_H

class BinMatrix {
  public:
    BinMatrix(int ydim,int xdim)
    :data(ydim*xdim),m_ydim(ydim),m_xdim(xdim),m_ncells(xdim*ydim)
    {

    }
    int GetYDim(){return m_ydim;};
    int GetXDim(){return m_xdim;};
    int GetNCells(){return m_ncells;};
    void Set(int y,int x,int val)
    {
      data[y*m_xdim+x]=val;
    }
    int Get(int y,int x)
    {
      return data[y*m_xdim+x];
    }
    void Print()
    {
       for (int j=0;j<m_ydim;j++)
       {
         for (int i=0;i<m_xdim;i++) cout << std::setw(2) << data[j*m_xdim+i] << " ";
         cout << endl;
       }
    }
  protected:
    std::vector <int> data;
    int m_ydim,m_xdim,m_ncells;
};

// extremly simple cluster labeling algorithm
class ClusterLabel {
  public:
    ClusterLabel(BinMatrix &Matrix,int len)
    :m(Matrix),ncells(m.GetXDim()*m.GetYDim()),total_edge_len(0),clusters((ncells+1)/2),m_len((double)len)
    {};
    void LabelClusters()
    {
        int max_cluster_label=0;
        for (int j=0;j<m.GetYDim();j++)
          for (int i=0;i<m.GetXDim();i++)
          {
             if (m.Get(j,i)) { // site is occupied
                int left,top,right,bottom;
                left=top=right=bottom=0;
                if (i>0) left=m.Get(j,i-1);
                if (j>0) top=m.Get(j-1,i);
                if (i<m.GetXDim()-1) right=m.Get(j,i+1);
                if (j<m.GetYDim()-1) bottom=m.Get(j+1,i);

                if (!left) total_edge_len++;
                if (!right) total_edge_len++;
                if (!top) total_edge_len++;
                if (!bottom) total_edge_len++;

                if (left==0 && top==0) {
                    max_cluster_label++;
                    m.Set(j,i,max_cluster_label);
                } else if (top==0 || left==0 || (top==left)) //cluster already seen
                {
                  m.Set(j,i,std::max(left,top));
                } else // conflict: map the larger cluster to the smaller one
                {
                  int cmin=std::min(left,top);
                  int cmax=std::max(left,top);
                  m.Set(j,i,cmin);
                  if (cmin!=cmax) {
                    for (int k=0;k<j;k++) RelabelRow(k,m.GetXDim(),cmin,cmax);
                    RelabelRow(j,i,cmin,cmax);
                  }
                }
             }
          }
      CalculateStats();
    }
    void RelabelRow(int row,int max_col,int cmin,int cmax)
    {
       for (int i=0;i<max_col;i++)
        if (m.Get(row,i)==cmax) m.Set(row,i,cmin);
    }
    cluster_stats &GetClusterStats(){return myStats;};
  protected:
    void CalculateStats()
    {
        for (int j=0;j<m.GetYDim();j++)
          for (int i=0;i<m.GetXDim();i++)
          {
            int label=m.Get(j,i);
            if (label)
            {
              clusters[label]++;
              myStats.cell_area++;
              myStats.total_area+=(m_len*m_len);
            }
          }
        for (size_t i=0;i<clusters.size();i++)
        {
           if (clusters[i])
           {
              double area=(double)clusters[i]*(m_len*m_len);
              myStats.num_clusters++;
             if (area>myStats.max_area) myStats.max_area=area;
           }
        }
        myStats.total_border_len=(double)total_edge_len*m_len;
    }
    cluster_stats myStats;
    BinMatrix &m;
    int ncells,total_edge_len;
    std::vector <int> clusters;
    double m_len;
};

#endif // BMATRIX_H
