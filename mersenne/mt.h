#define kN 624   /* Period parameters */
#define kM 397

extern inline unsigned long int mt_get (void *vstate);
extern long double mt_get_double (void *vstate);
extern void mt_set (void *state, unsigned long int s);

typedef struct
  {
    unsigned long mt[kN];
    int mti;
  }mt_state_t;

#undef kN
#undef kM
