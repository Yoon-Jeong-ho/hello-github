#include <iostream>
#include<string>


using namespace std;

int main()
{
    string str;
    int result = 1;
    getline(cin,str);
    if(str.length()==1 && str[0]==' '){
        cout<<0;
        return 0;
    }

    for(int i = 1 ; i < str.length() ; i++){
        if(str[i]==' ')result++;
    }
    cout<<result;
    return 0;
}
