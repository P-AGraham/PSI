

template<typename Indtype, typename Dtype>
void quick_sort(Indtype* indices, Dtype* data, size_t l, size_t r)// sort data by indices, r = length-1
{
    if (l < r)
    {
        size_t i,j;
        Indtype x;
        Dtype y;

        i = l;
        j = r;
        x = indices[i];
        y = data[i];
        while (i < j)
        {
            while(i < j && indices[j] > x)
            {
                j--; // 从右向左找第一个小于x的数
            }
            if(i < j)
            {
                indices[i] = indices[j];
                data[i] = data[j];
                i++;
            }
            while(i < j && indices[i] < x)
            {
                i++; // 从左向右找第一个大于x的数
            }
            if(i < j)
            {
                indices[j] = indices[i];
                data[j] = data[i];
                j--;
            }
        }
        indices[i] = x;
        data[i] = y;
        if(i>l)
            quick_sort(indices, data, l, i-1); /* 递归调用 */
        if(r>i)
            quick_sort(indices, data, i+1, r); /* 递归调用 */
    }
}