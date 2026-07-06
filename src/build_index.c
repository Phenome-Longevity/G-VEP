/* build_index: read uncompressed cache stream on stdin + the bgzip .gzi,
 * emit one 16-byte record per line: (FNV1a-64 of the key) + (bgzf virtual offset).
 * Usage: bgzip -dc cache.bgz | ./build_index cache.bgz.gzi > index_raw.bin
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(int argc, char **argv){
    if(argc<2){ fprintf(stderr,"usage: build_index <gzi>\n"); return 1; }
    FILE *g = fopen(argv[1],"rb");
    if(!g){ perror("gzi"); return 1; }
    uint64_t n;
    if(fread(&n,8,1,g)!=1){ fprintf(stderr,"bad gzi\n"); return 1; }
    /* blocks 0..n ; block 0 = (0,0), block i = gzi entry i-1 */
    uint64_t *co = malloc((n+1)*8), *uo = malloc((n+1)*8);
    co[0]=0; uo[0]=0;
    for(uint64_t i=1;i<=n;i++){ if(fread(&co[i],8,1,g)!=1||fread(&uo[i],8,1,g)!=1){fprintf(stderr,"gzi trunc\n");return 1;} }
    fclose(g);

    char *line=NULL; size_t cap=0; ssize_t len;
    uint64_t U=0, b=0, count=0;
    /* big output buffer */
    setvbuf(stdout, NULL, _IOFBF, 1<<22);
    while((len=getline(&line,&cap,stdin))>=0){
        while(b+1<=n && uo[b+1]<=U) b++;
        uint64_t voff = (co[b]<<16) | (U - uo[b]);
        char *tab = memchr(line,'\t',len);
        size_t klen = tab ? (size_t)(tab-line) : (size_t)len;
        uint64_t h = 14695981039346656037ULL;
        for(size_t i=0;i<klen;i++){ h ^= (unsigned char)line[i]; h *= 1099511628211ULL; }
        fwrite(&h,8,1,stdout); fwrite(&voff,8,1,stdout);
        U += (uint64_t)len;
        count++;
    }
    fprintf(stderr,"build_index: %llu records, %llu uncompressed bytes\n",
            (unsigned long long)count,(unsigned long long)U);
    return 0;
}
