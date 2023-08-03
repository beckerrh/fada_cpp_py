#include  "matrixinterface.hpp"

/*-------------------------------------------------*/
SystemMatrix::SystemMatrix(std::map<std::string,std::shared_ptr<MatrixInterface>> matrices, std::map<std::pair<int,int>,std::pair<std::string, std::string>> patterns):_matrices(matrices), _patterns(patterns) 
{
    int bsize=0;
    for(auto p:patterns)
    {
        bsize = std::max(bsize, p.first.first);
        bsize = std::max(bsize, p.first.second);
    }
    bsize += 1;
    // std::cerr << "bsize="<<bsize<<"\n";
    armaicvec n(bsize);
    n.fill(-1);
    for(auto p:patterns)
    {
        // std::cerr << "n="<<n.t();
        assert(matrices.at(p.second.first));
        if(p.second.second=="T" or p.second.second=="MT")
        {
            int i(p.first.first), j(p.first.second), nr(matrices.at(p.second.first)->ncols()), nc(matrices.at(p.second.first)->nrows());
            // std::cerr << "TT i="<<i << " j="<< j << " nr="<<nr << " nc=" << nc << "\n";
            if(n[i]!=-1)
            {
                assert(n[i]==nr);
            }
            else
            {
                n[i] = nr;
            }
            if(n[j]!=-1)
            {
                assert(n[j]==nc);
            }
            else
            {
                n[j] = nc;
            }
            // n[p.first.first] = matrices.at(p.second.first)->ncols();
            // n[p.first.second] = matrices.at(p.second.first)->nrows();
        }
        else
        {
            int i(p.first.first), j(p.first.second), nc(matrices.at(p.second.first)->ncols()), nr(matrices.at(p.second.first)->nrows());
            // std::cerr << "-- i="<<i << " j="<< j << " nr="<<nr << " nc=" << nc << "\n";
            if(n[i]!=-1)
            {
                assert(n[i]==nr);
            }
            else
            {
                n[i] = nr;
            }
            if(n[j]!=-1)
            {
                assert(n[j]==nc);
            }
            else
            {
                n[j] = nc;
            }
            // n[p.first.first] = matrices.at(p.second.first)->nrows();
            // n[p.first.second] = matrices.at(p.second.first)->ncols();
        }
    }
    // std::cerr << "n="<<n<<"\n";
    _ofs.set_size(bsize+1);
    _ofs[0] = 0;
    for(int i=0;i<bsize;i++)
    {
        _ofs[i+1] = _ofs[i] + n[i];
    }
    // std::cerr << "_ofs="<<n<<"\n";
}
