/* lookup: stdin = sample keys (one per line). Loads sorted (hash->voffset) index (mmap)
 * + cache.bgz. For each key: FNV1a-64 -> binary search -> bgzf seek+read -> verify key ->
 * print the cache line (key\tsuffix). Misses print nothing (they go to live-VEP fallback).
 * Usage: ./lookup h_sorted.bin v_sorted.bin cache.bgz < keys > sample_cache.tsv
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

static inline uint64_t fnv(const char *s, size_t n){
    uint64_t h=14695981039346656037ULL;
    for(size_t i=0;i<n;i++){ h ^= (unsigned char)s[i]; h *= 1099511628211ULL; }
    return h;
}
static uint64_t* mapu64(const char *path, size_t *count){
    int fd=open(path,O_RDONLY);
    if(fd<0){ perror(path); exit(1);}
    struct stat st; fstat(fd,&st);
    uint64_t *p=mmap(NULL,st.st_size,PROT_READ,MAP_SHARED,fd,0);
    if(p==MAP_FAILED){ perror("mmap"); exit(1);}
    close(fd); *count=st.st_size/8; return p;
}
int main(int argc,char**argv){
    if(argc<4){ fprintf(stderr,"usage: lookup h_sorted v_sorted cache.bgz\n"); return 1; }
    size_t N,Nv; uint64_t *H=mapu64(argv[1],&N); uint64_t *V=mapu64(argv[2],&Nv);
    BGZF *fp=bgzf_open(argv[3],"r");
    if(!fp){ fprintf(stderr,"cannot open %s\n",argv[3]); return 1; }
    kstring_t ks={0,0,0};
    char *line=NULL; size_t cap=0; ssize_t len;
    uint64_t hits=0,total=0;
    setvbuf(stdout,NULL,_IOFBF,1<<22);
    while((len=getline(&line,&cap,stdin))>=0){
        total++;
        while(len>0 && (line[len-1]=='\n'||line[len-1]=='\r')) line[--len]=0;
        if(len==0) continue;
        uint64_t h=fnv(line,(size_t)len);
        size_t lo=0,hi=N;
        while(lo<hi){ size_t mid=lo+((hi-lo)>>1); if(H[mid]<h) lo=mid+1; else hi=mid; }
        for(size_t i=lo;i<N && H[i]==h;i++){
            if(bgzf_seek(fp,(int64_t)V[i],SEEK_SET)<0) continue;
            ks.l=0;
            if(bgzf_getline(fp,'\n',&ks)<0) continue;
            char *tab=memchr(ks.s,'\t',ks.l);
            size_t klen=tab?(size_t)(tab-ks.s):ks.l;
            if(klen==(size_t)len && memcmp(ks.s,line,(size_t)len)==0){
                fwrite(ks.s,1,ks.l,stdout); fputc('\n',stdout);
                hits++; break;
            }
        }
    }
    fprintf(stderr,"lookup: %llu/%llu hits\n",(unsigned long long)hits,(unsigned long long)total);
    return 0;
}
