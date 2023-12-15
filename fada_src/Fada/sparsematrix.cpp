//
//  sparsematrix.cpp
//  Fada
//
//  Created by Roland Becker on 20/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <set>
#include  "sparsematrix.hpp"
#include  "matrixinterface.hpp"
#include  "systemmatrix.hpp"

/*-------------------------------------------------*/
void SparseMatrix::save(std::ostream& out, arma::file_type datatype) const
{
    for(int i=0;i<nrows();i++)
    {
        for(int pos=_rows[i];pos<_rows[i+1];pos++)
        {
            out << "(" << i << " " << _cols[pos] << ") " << _values[pos] << "\n";
        }
    }
}
/*-------------------------------------------------*/
void SparseMatrix::dot(armavec& x, const armavec& b, double d) const
{
    assert(nrows()==x.n_elem);
    assert(ncols()==b.n_elem);
    for(int i=0;i<nrows();i++)
    {
        for(int pos=_rows[i];pos<_rows[i+1];pos++)
        {
            x[i] += d*_values[pos] * b[_cols[pos]];
        }
    }
}
/*-------------------------------------------------*/
void SparseMatrix::Tdot(armavec& x, const armavec& b, double d) const
{
    assert(nrows()==b.n_elem);
    assert(ncols()==x.n_elem);
    for(int i=0;i<nrows();i++)
    {
        for(int pos=_rows[i];pos<_rows[i+1];pos++)
        {
            x[_cols[pos]] += d*_values[pos] * b[i];
        }
    }
}
/*-------------------------------------------------*/
void SparseMatrix::row_zero(int i)
{
    for(int pos=_rows[i];pos<_rows[i+1];pos++)
    {
        _values[pos]=0.0;
    }    
}

/*-------------------------------------------------*/
std::shared_ptr<armavec> SparseMatrix::getDiag() const
{
    auto diag = std::make_shared<armavec>(nrows());
    // armavec diag(nrows());
    for(int i=0;i<nrows();i++)
    {
        // diag[i] = _values[_diags[i]];
        diag->at(i) = _values[_diags[i]];
    }
    return diag;
}
/*-------------------------------------------------*/
std::shared_ptr<SparseMatrix> SparseMatrix::getTransposed(double d) const
{
    arma::umat locations(2,_values.n_elem);
    int count=0;
    for(int i=0;i<nrows();i++)
    {
        for(int pos=_rows[i];pos<_rows[i+1];pos++)
        {
            locations.at(0,count) = _cols[pos];
            locations.at(1,count) = i;
            count++;
        }
    }    
    return std::make_shared<SparseMatrix>(locations, d*_values, _diags.n_elem==nrows());
}
/*-------------------------------------------------*/
void SparseMatrix::from_transposed(const SparseMatrix& matrix, double d)
{
    arma::umat locations(2,matrix.nelem());
    int count=0;
    for(int i=0;i<matrix.nrows();i++)
    {
        for(int pos=matrix.rows()[i];pos<matrix.rows()[i+1];pos++)
        {
            locations.at(0,count) = matrix.cols()[pos];
            locations.at(1,count) = i;
            count++;
        }
    }
    set_elements(locations, d*matrix.values(), matrix.diags().n_elem);        
}


/*-------------------------------------------------*/
void SparseMatrix::compute_diag()
{
    _diags.set_size(nrows());
    for(int i=0;i<nrows();i++)
    {
        // bool found=false;
        for(int pos=_rows[i];pos<_rows[i+1];pos++)
        {
            if(_cols[pos]==i)
            {
                _diags[i] = pos;
                // found = true;
                break;
            }
        }
        // if (not found)
        // {
        //     for(int pos=_rows[i];pos<_rows[i+1];pos++)
        //     {
        //         std::cerr << "i="<<i << "pos="<<pos << "_cols[pos]="<<_cols[pos]  << "_rows[i]="<<_rows[i] << "_rows[i+1]=" << _rows[i+1] <<"\n";
        //     }
        //     assert(0);
        // }
    }
    // for(int i=0;i<nrows;i++)
    // {
    //   std::cerr << "\nrow = " << i << " " << _cols[_diags[i]] << "\n";
    //   for(int pos=_rows[i];pos<_rows[i+1];pos++)
    //   {
    //     std::cerr << " " << _cols[pos];
    //   }
    // }
    // save(std::cerr);
    //  assert(0);    
}
/*-------------------------------------------------*/
void SparseMatrix::add_diagonal(double d)
{
    assert(_diags.n_rows==nrows());
    for(int i=0;i<nrows();i++)
    {
        _values[_diags[i]] += d;
    }
}


/*-------------------------------------------------*/
void SparseMatrix::set_elements(const arma::umat& locations, const armavec& values, bool compute_diag_flag)
{
    // No duplicates in locations !!!!
    assert(locations.n_rows==2);
    assert(locations.n_cols==values.n_elem);
    // _values.fill(0);
    // std::cerr << "locations\n" << locations.row(0) << "\n" << locations.row(1) << "\n";
    // sort the indices !!
    _n_cols = arma::max(locations.row(1))+1;
    arma::uvec ind = arma::sort_index(_n_cols*locations.row(0)+locations.row(1));
    // CHECK FOR REPTITIONS
    // arma::uvec un = arma::unique(ind);
    // assert(un.n_elem==_n_cols);    
    arma::uword n = locations.n_cols;
    arma::uword nrows = locations(0,ind[n-1])+1;
    // std::cerr << "nrows=" << nrows << "\n";
    // std::cerr << "_n_cols=" << _n_cols << "\n";
    // std::cerr << "arma::max(locations.row(1))*locations.row(0)+locations.row(1)=" << arma::max(locations.row(1))*locations.row(0)+locations.row(1) << "\n";
    // std::cerr << "ind=" << ind.t() << "\n";
    // for(int i=0;i<n;i++)
    // {
    //   arma::uword index = ind[i];
    //   arma::uword row = locations.at(0,index);
    //   std::cerr << "i="<<i << " index="<<index << " row=" <<row << " col=" <<locations.at(1,index) << "\n";
    //   }
    bool allow_duplicates=false;
    if(allow_duplicates)
    {
        arma::uvec un = arma::unique(ind);
        arma::uword nreal = un.n_elem;
        _cols.set_size(nreal);
        _values.set_size(nreal);
        _values.fill(0);
        _rows.set_size(nrows+1);
        int count = 0, prec=-1;
        int count_col = 0, prec_col=-1;
        for(int i=0;i<n;i++)
        {
            arma::uword index = ind[i];
            arma::uword row = locations.at(0,index);
            // std::cerr << "i="<<i << " row="<<row << " prec=" << prec << "\n";
            if(row!=prec)
            {
                _rows[count] = i;
                prec = row;
                count++;
                prec_col=-1;
            }
            arma::uword col = locations.at(1,index);
            if(col!=prec_col)
            {
                _cols[count_col] = col;
                _values[count_col] += values[index];
                prec_col = col;
                count_col++;
            }
            else
            {
                _values[count_col] += values[index];                
            }
        }
        _rows[nrows] = nreal;
        
    }
    else
    {
        _cols.set_size(n);
        _values.set_size(n);
        _rows.set_size(nrows+1);
        int count = 0, prec=-1;
        for(int i=0;i<n;i++)
        {
            arma::uword index = ind[i];
            arma::uword row = locations.at(0,index);
            // std::cerr << "i="<<i << " row="<<row << " prec=" << prec << "\n";
            if(row!=prec)
            {
                _rows[count] = i;
                prec = row;
                count++;
            }
            _values[i] = values[index];
            _cols[i] = locations.at(1,index);
        }
        _rows[nrows] = n;        
    }
    // std::cerr << "_rows=" << _rows.t() << "\n";
    if(compute_diag_flag)
    {
        compute_diag();
    }
}
/*-------------------------------------------------*/
void SparseMatrix::jacobi(armavec& out, const armavec& in) const
{
    for(int i=0;i<nrows();i++)
    {
        out[i] = in[i]/_values[_diags[i]];
    }
}
/*-------------------------------------------------*/
void SparseMatrix::gauss_seidel1(armavec& out, const armavec& in) const
{
    int niter=2;
    // std::cerr << "gauss_seidel1 " << "in=" << in.t()<<"\n";
    for(int iter=0;iter<niter;iter++)
    {     
        for(int i=0;i<nrows();i++)
        {
            double val = in[i];
            for(int pos=_rows[i];pos<_rows[i+1];pos++)
            {
                int j = _cols[pos];
                // if(j>=i) break;
                if(i!=j)
                    val -= _values[pos] * out[j];
            }
            out[i] = val/_values[_diags[i]];
        }
    }
    // std::cerr << "gauss_seidel1 " << "out=" << out.t()<<"\n";
}
/*-------------------------------------------------*/
void SparseMatrix::gauss_seidel2(armavec& out, const armavec& in) const
{
    int niter=2;
    // std::cerr << "gauss_seidel2 " << "nrows()=" << nrows() << " in.n_elem=" << in.n_elem<<"\n";
    // std::cerr << "gauss_seidel2 " << "in=" << in.t()<<"\n";
    // std::cerr << "gauss_seidel2 " << "diag=" << getDiag().t()<<"\n";
    for(int iter=0;iter<niter;iter++)
    {
        for(int i=nrows()-1;i>=0;i--)
        {
            double val = in[i];
            for(int pos=_rows[i+1]-1;pos>=_rows[i];pos--)
            {
                int j = _cols[pos];
                // if(j<=i) break;
                if(i!=j)
                    val -= _values[pos] * out[j];
            }
            out[i] = val/_values[_diags[i]];
        }
    }
    // std::cerr << "gauss_seidel2 " << "out=" << out.t()<<"\n";
}
/*-------------------------------------------------*/
void SparseMatrix::addBBT(const SparseMatrix& B, double d, std::shared_ptr<armavec const> D)
{
    // std::cerr << "B.nrows() " << B.nrows() << " B.ncols() " << B.ncols() <<"\n";

#ifdef _LONG_LONG
    const arma::uvec& Brows = B.rows();
    const arma::uvec& Bcols = B.cols();
#else
    const armaicvec& Brows = B.rows();
    const armaicvec& Bcols = B.cols();
#endif
    const armavec& Bvalues = B.values();

    if(D)
    {
        assert(D->n_elem==B.ncols());
    }
    std::vector<std::set<int>> rows_inv(B.ncols());
    for(int i=0;i<B.nrows();i++)
    {
        for(int posi=Brows[i];posi<Brows[i+1];posi++)
        {
            int k = Bcols[posi];
            rows_inv[k].insert(i);
        }
    }
        
    std::vector<std::set<int>> rows(B.nrows());
    // insert old values
    for(int i=0;i<nrows();i++)
    {
        for(int pos=_rows[i];pos<_rows[i+1];pos++)
        {
            rows[i].insert(_cols[pos]);
        }            
    }
    for(int i=0;i<B.nrows();i++)
    {
        for(int posi=Brows[i];posi<Brows[i+1];posi++)
        {
            int k = Bcols[posi];
            rows[i].insert(rows_inv[k].begin(), rows_inv[k].end());
        }
    }
        
#ifdef _LONG_LONG
    arma::uvec rows_save(_rows);
    arma::uvec cols_save(_cols);
#else
    armaicvec rows_save(_rows);
    armaicvec cols_save(_cols);
#endif
    armavec values_save(_values);
        
        
    _rows.resize(B.nrows()+1);
    int n=0;
    _rows[0] = 0;
    for(int i=0;i<rows.size();i++)
    {
        int ni = rows[i].size();
        _rows[i+1] = _rows[i] + ni;
        n += ni;
    }
    _cols.resize(n);
    int count=0;
    for(int i=0;i<rows.size();i++)
    {
        for(auto p:rows[i])
        {
            _cols[count++] = p;
        }
    }
    _values.resize(n);
    _values.fill(0.0);
        
    
    // copy old values
    for(int i=0;i<rows_save.n_elem-1;i++)
    {
        for(int pos=rows_save[i];pos<rows_save[i+1];pos++)
        {
            int j = cols_save[pos];
            for(int pos2=_rows[i];pos2<_rows[i+1];pos2++)
            {
                if(_cols[pos2]==j)
                {
                    _values[pos2] += values_save(pos);
                }
            }            
        }
    }            
    for(int i=0;i<B.nrows();i++)
    {
        double scale=d;
        if(D) scale /= D->at(i);
        for(int pos=_rows[i];pos<_rows[i+1];pos++)
        {
            int j = _cols[pos];
            for(int posi=Brows[i];posi<Brows[i+1];posi++)
            {
                int k = Bcols[posi];
                for(int posj=Brows[j];posj<Brows[j+1];posj++)
                {
                    int k2 = Bcols[posj];
                    if(k2==k)
                    {
                        // posj = (j,k)
                        _values[pos] += scale*Bvalues[posi]*Bvalues[posj];
                    }
                }   
            }
        }
    }
    compute_diag();
    // B.save(std::cerr);
    // std::cerr << "and now:\n";
    // this->save(std::cerr);
}

/*-------------------------------------------------*/
void SparseMatrix::set_elements(std::shared_ptr<SystemMatrix const> sm)
{
    const std::map<std::string,std::shared_ptr<MatrixInterface>>& matrices =  sm->get_matrices();
    const std::map<std::pair<int,int>,std::pair<std::string, std::string>>& patterns = sm->get_patterns();
    const armaicvec& ofs = sm->get_ofs();
    int n = ofs.n_elem-1;
    // std::cerr << "n="<<n<<"\n";
    std::vector<std::vector<std::string>> mats(n,std::vector<std::string>(n,"")), types(n,std::vector<std::string>(n,""));
    for(auto p:patterns)
    {
        mats[p.first.first][p.first.second] = p.second.first;
        types[p.first.first][p.first.second] = p.second.second;
    }
    // for(int i=0;i<n;i++)
    // {
    //     for(int j=0;j<n;j++)
    //     {
    //         std::cerr << "i j " << i << " " << j << " " << mats[i][j] << " " << types[i][j] << "\n";
    //     }
    // }
    int count=0;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(mats[i][j]!="")
            {
                auto sp = std::dynamic_pointer_cast<SparseMatrix const>(matrices.at(mats[i][j]));
                int nextsize = sp->nelem();
                // std::cerr << "nextsize="<<nextsize<<" " << mats[i][j] << "\n";
                count += nextsize;
            }
        }
    }
    // std::cerr << "count="<<count<<"\n";
    armavec values(count);
    arma::umat locations(2,count);
    count=0;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(mats[i][j]!="")
            {
                auto sp = std::dynamic_pointer_cast<SparseMatrix const>(matrices.at(mats[i][j]));
                int nextsize = sp->nelem();
                // std::cerr << "nextsize="<<nextsize<<" " << mats[i][j] << " " << count << "\n";
                if(types[i][j]=="MT")
                {
                    values(arma::span(count,count+nextsize-1)) = -sp->values();
                }
                else
                {
                    values(arma::span(count,count+nextsize-1)) = sp->values();                    
                }
                count += nextsize;
            }
        }
    }
    count=0;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            // std::cerr << "i j " << i << " " << j << " " << mats[i][j] << " " << types[i][j] << "\n";
            if(mats[i][j]!="")
            {
                auto sp = std::dynamic_pointer_cast<SparseMatrix const>(matrices.at(mats[i][j]));
                for(int ii=0;ii<sp->nrows();ii++)
                {
                    for(int pos=sp->rows()[ii];pos<sp->rows()[ii+1];pos++)
                    {
                        int jj = sp->cols()[pos];
                        if(types[i][j]=="T" or types[i][j]=="MT")
                        {
                            // std::cerr << "i="<<i << " j="<<j <<  " ofs[i]="<<ofs[i] << " ofs[j]="<<ofs[j] << " ii="<<ii<<" jj="<<jj<<"\n";
                            locations.at(0,count) = jj + ofs[i];
                            locations.at(1,count) = ii + ofs[j]; 
                        }
                        else
                        {
                            locations.at(0,count) = ii + ofs[i];                             
                            locations.at(1,count) = jj + ofs[j]; 
                        }
                        count++; 
                    }
                }
                // std::cerr << "locations\n" << locations.t();
            }
        }
    }
    set_elements(locations, values);    
}
