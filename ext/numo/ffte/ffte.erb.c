/*
  FFTE: A Fast Fourier Transform Package
  http://www.ffte.jp/
  Computes Discrete Fourier Transforms of
  1-, 2- and 3- dimensional sequences of length (2^p)*(3^q)*(5^r).
  Copyright(C) 2000-2004,2008-2011 Daisuke Takahashi

  Ruby/NArray wrapper of FFTE
  FFTE is converted to C-code using f2c.
  Copyright(C) 2013 Masahiro Tanaka
*/

#include <ruby.h>
#include "numo/narray.h"
#include "numo/template.h"
//#include "f2c.h"
#define cDC numo_cDComplex
#define cDF numo_cDFloat

typedef int integer;

static VALUE eRadixError;

static inline
int is235radix(integer n) {
    if (n<=1) {return 0;}
    while (n % 5 == 0) {n /= 5;}
    while (n % 3 == 0) {n /= 3;}
    return (n & (n-1)) ? 0 : 1;
}

typedef struct {
    dcomplex *buf;
    integer iopt;
    integer dummy;
    size_t nb;
} fft_opt_t;

static void
fft_opt_free(void *ptr)
{
    fft_opt_t *g = (fft_opt_t*)ptr;
    xfree(g->buf);
    xfree(g);
}

static size_t
fft_opt_memsize(const void *ptr)
{
    const fft_opt_t *g = (const fft_opt_t*)ptr;
    return sizeof(fft_opt_t) + g->nb * sizeof(dcomplex);
}

static const rb_data_type_t fft_opt_type = {
    "Numo::FFTE/opt",
    {NULL, fft_opt_free, fft_opt_memsize,},
    0, 0, RUBY_TYPED_FREE_IMMEDIATELY|RUBY_TYPED_WB_PROTECTED
};

static inline fft_opt_t *
alloc_fft_opt(int nb, integer iopt, VALUE *v)
{
    fft_opt_t *g;

    g = ALLOC(fft_opt_t);
    g->buf = ALLOC_N(dcomplex, nb);
    g->iopt = iopt;
    g->nb = nb;
    *v = TypedData_Wrap_Struct(rb_cData, &fft_opt_type, (void*)g);
    return g;
}

<%
def argmap(arg)
  case arg
  when Numeric
    arg = 1..arg
  end
  arg.map do |x|
    yield(x)
  end.join(",")
end
$funcs = []
%>
int zfft1d_(dcomplex *a, integer *n, integer *iopt, dcomplex *b);
int zfft2d_(dcomplex *a, integer *nx, integer *ny, integer *iopt);
int zfft3d_(dcomplex *a, integer *nx, integer *ny, integer *nz, integer *iopt);

<% (1..3).each do |d| %>
// <%=d%>-dimentional 2,3,4,5,8-radix FFT
static void
iter_fft_zfft<%=d%>d(na_loop_t *const lp)
{
    char   *p1;
    integer iopt;
    integer <%=argmap(d){|i|"n#{i}"}%>;
<% if d==1 %>
    fft_opt_t *g;
<% end %>

<% (1..d).each do |i| %>
    //n<%=i%> = lp->n[<%=i-1%>];
    n<%=i%> = lp->args[0].shape[<%=i-1%>];
<% end %>
    p1 = ((lp)->args[0]).ptr + ((lp)->args[0].iter[0]).pos;

<% if d==1 %>
    g = (fft_opt_t*)(lp->opt_ptr);
    iopt = g->iopt;
    zfft<%=d%>d_((dcomplex*)p1, <%=argmap(d){|i|"&n#{i}"}%>, &iopt, g->buf);
<% else %>
    iopt = *(integer*)(lp->opt_ptr);
    zfft<%=d%>d_((dcomplex*)p1, <%=argmap(d){|i|"&n#{i}"}%>, &iopt);
<% end %>
}


/*
  <%=d%>-dimentional COMPLEX FFT
  using Radix-2,3,4,5,8 FFT routine.
  Calculates on each last <%=d%>-dimention.
  @overload zfft<%=d%>d(narray,[iopt])
  @param [Numo::DComplex] narray \>=<%=d%>-dimentional Complex NArray with
    `shape = [.., <%if d>2%>nz, <%end;if d>1%>ny, <%end%>nx]` where
    `nx = (2**ip) * (3**iq) * (5**ir)<% if d>1 %>`,
    `ny = (2**jp) * (3**jq) * (5**jr)<% if d>2 %>`,
    `nz = (2**kp) * (3**kq) * (5**kr)<% end; end %>`.
  @param [Numeric] iopt  Transform direction.
    `iopt=-1` for FORWARD transform,
    `iopt=+1` for INVERSE transform.
    Default: `iopt=+1`.
  @return [Numo::DComplex] Result DComplex NArray with
    `shape = [.., <%if d>2%>nz, <%end;if d>1%>ny, <%end%>nx]`.
  @raise  [Numo::FFTE::RadixError] if `nx`<%if d>1%>, `ny`<% if d>2%>, `nz`<%end;end%>
    is not `(2^p)*(3^q)*(5^r)`.
*/
<% $funcs.push func="zfft#{d}d" %>
static VALUE
numo_ffte_<%=func%>(int argc, VALUE *args, VALUE mod)
{
    narray_t *na;
    VALUE vres, viopt;
    VALUE vna;
    int ndim;
    integer iopt=1, iopt_zero=0;
    ndfunc_arg_in_t ain[1] = {{cDC,<%=d%>}};
    ndfunc_t ndf = { iter_fft_zfft<%=d%>d, NO_LOOP, 1, 0, ain, 0 };
<% if d==1 %>
    fft_opt_t *g;
    VALUE vopt;
<% end %>
    integer <%=argmap(d){|i|"n#{i}"}%>;

    switch(rb_scan_args(argc, args, "11", &vna, &viopt)) {
    case 2:
        iopt = NUM2INT(viopt);
    }
    GetNArray(vna,na);
    ndim = NA_NDIM(na);
    if (ndim<<%=d%>) {
        rb_raise(eRadixError,"ndim(=%d) should >= <%=d%>",ndim);
    }
<% (1..d).each do |i| %>
    n<%=i%> = NA_SHAPE(na)[NA_NDIM(na)-<%=i%>];
    if (!is235radix(n<%=i%>)) {
        rb_raise(eRadixError,"%d-th dim length(=%d) is not 2,3,5-radix",NA_NDIM(na)-<%=i%>,n<%=i%>);
    }
<% end %>

    vres = na_copy(vna);

<% if d==1 %>
    g = alloc_fft_opt(n1*2, iopt, &vopt);
    zfft1d_(NULL, &n1, &iopt_zero, g->buf);
    na_ndloop3(&ndf, g, 1, vres);
    RB_GC_GUARD(vopt);
<% else %>
    zfft<%=d%>d_(NULL,  <%=argmap(d){|i|"&n#{i}"}%>, &iopt_zero);
    na_ndloop3(&ndf, &iopt, 1, vres);
<% end %>

    return vres;
}
<% end %>


int zdfft2d_(dcomplex *a, integer *nx, integer *ny, integer *iopt, dcomplex *b);
int zdfft3d_(dcomplex *a, integer *nx, integer *ny, integer *nz, integer *iopt, dcomplex *b);

<% (2..3).each do |d| %>
// <%=d%>-dimentional 2,3,4,5,8-radix FFT
static void
iter_fft_zdfft<%=d%>d(na_loop_t *const lp)
{
    char *p1, *p2;
    integer <%=argmap(d){|i|"n#{i}"}%>;
    size_t i, n;
    fft_opt_t *g;

    //n1 = n = (lp->n[<%=d-1%>]-1)*2;
    n1 = n = (lp->args[0].shape[<%=d-1%>]-1)*2;
<% (2..d).each do |i| %>
    //n<%=i%> = lp->n[<%=d-i%>];
    n<%=i%> = lp->args[0].shape[<%=d-i%>];
    n *= n<%=i%>;
<% end %>
    g = (fft_opt_t*)(lp->opt_ptr);
    p1 = ((lp)->args[0]).ptr + ((lp)->args[0].iter[0]).pos;
    p2 = ((lp)->args[1]).ptr + ((lp)->args[1].iter[0]).pos;

    zdfft<%=d%>d_((dcomplex*)p1, <%=argmap(d){|i|"&n#{i}"}%>, &(g->iopt), g->buf);

    for (i=0; i<n; i++) {
        ((double*)p2)[i] = ((double*)p1)[i];
    }
}


/*
  <%=d%>-dimentional COMPLEX-TO-REAL FFT
  using Radix-2,3,4,5,8 FFT routine.
  Calculates on each last <%=d%>-dimention.
  This routine is INVERSE transform.
  @overload zdfft<%=d%>d(narray)
  @param [Numo::DComplex] narray
    \>=<%=d%>-dimentional DComplex NArray with
    `shape =[.., <%if d>2%>nz, <%end;if d>1%>ny, <%end%>nx/2+1]` where
    `nx = (2**ip) * (3**iq) * (5**ir)`,
    `ny = (2**jp) * (3**jq) * (5**jr)<% if d==3 %>`,
    `nz = (2**kp) * (3**kq) * (5**kr)<% end %>`.
  @return [Numo::DFloat]  Result DFloat NArray with
    `shape = [.., <%if d>2%>nz, <%end;if d>1%>ny, <%end%>nx]`
  @raise  [Numo::FFTE::RadixError] if `nx`, `ny`<%if d>2%>, `nz`<%end%>
    is not `(2^p)*(3^q)*(5^r)`.
*/
<% $funcs.push func="zdfft#{d}d" %>
static VALUE
numo_ffte_<%=func%>(int argc, VALUE *args, VALUE mod)
{
    narray_t *na;
    VALUE vres;
    VALUE vb, vna, viopt;
    int ndim;
    integer iopt=1, iopt_zero=0;
    fft_opt_t *g;
    integer <%=argmap(d){|i|"n#{i}"}%>;
    size_t n=1;
    size_t shape[<%=d%>];
    ndfunc_arg_in_t ain[1] = {{cDC,<%=d%>}};
    ndfunc_arg_out_t aout[1] = {{cDF,<%=d%>,shape}};
    ndfunc_t ndf = { iter_fft_zdfft<%=d%>d, NO_LOOP, 1, 1, ain, aout };

    switch(rb_scan_args(argc, args, "11", &vna, &viopt)) {
    case 2:
        iopt = NUM2INT(viopt);
    }
    GetNArray(vna,na);
    ndim = NA_NDIM(na);
    if (ndim<<%=d%>) {
        rb_raise(eRadixError,"ndim(=%d) should >= <%=d%>",ndim);
    }
<% (1..d).each do |i| %>
    n<%=i%> = NA_SHAPE(na)[NA_NDIM(na)-<%=i%>];
    n *= n<%=i%>;
<% if i==1 %>
    shape[<%=d-1%>] = (n<%=i%>-1)*2;
<% else %>
    shape[<%=d-i%>] =  n<%=i%>;
<% end %>
    if (!is235radix(shape[<%=d-i%>])) {
        rb_raise(eRadixError,"%d-th dim length(=%ld) is not 2,3,5-radix",
                 NA_NDIM(na)-<%=i%>,shape[<%=d-i%>]);
    }
<% end %>

    vna = na_copy(vna);
    g = alloc_fft_opt(n, iopt, &vb);

    zdfft<%=d%>d_(NULL, <%=argmap(d){|i|"&n#{i}"}%>, &iopt_zero, g->buf);
    vres = na_ndloop3(&ndf, g, 1, vna);
    RB_GC_GUARD(vb);

    return vres;
}
<% end %>


int dzfft2d_(dcomplex *a, integer *nx, integer *ny, integer *iopt, dcomplex *b);
int dzfft3d_(dcomplex *a, integer *nx, integer *ny, integer *nz, integer *iopt, dcomplex *b);

<% (2..3).each do |d| %>
// <%=d%>-dimentional 2,3,4,5,8-radix FFT
static void
iter_fft_dzfft<%=d%>d(na_loop_t *const lp)
{
    char *p1, *p2;
    integer <%=argmap(d){|i|"n#{i}"}%>;
    size_t i, n=1;
    fft_opt_t *g;

<% (1..d).each do |i| %>
    //n<%=i%> = lp->n[<%=d-i%>];
    n<%=i%> = lp->args[0].shape[<%=d-i%>];
    n *= n<%=i%>;
<% end %>
    g = (fft_opt_t*)(lp->opt_ptr);
    p1 = ((lp)->args[0]).ptr + ((lp)->args[0].iter[0]).pos;
    p2 = ((lp)->args[1]).ptr + ((lp)->args[1].iter[0]).pos;

    for (i=0; i<n; i++) {
        ((double*)p2)[i] = ((double*)p1)[i];
    }
    dzfft<%=d%>d_((dcomplex*)p2, <%=argmap(d){|i|"&n#{i}"}%>, &(g->iopt), g->buf);
}


/*
  <%=d%>-dimentional REAL-TO-COMPLEX FFT
  using Radix-2,3,4,5,8 FFT routine.
  Calculates on each last <%=d%>-dimention.
  This routine is FORWARD transform.
  @overload dzfft<%=d%>d(narray)
  @param [Numo::DFloat] narray \>=<%=d%>-dimentional DFloat NArray with
    `shape = [.., <%if d>2%>nz, <%end;if d>1%>ny, <%end%>nx]` where
    `nx = (2**ip) * (3**iq) * (5**ir)`,
    `ny = (2**jp) * (3**jq) * (5**jr)<% if d==3 %>`,
    `nz = (2**kp) * (3**kq) * (5**kr)<% end %>`.
  @return [Numo::DComplex]  Result DComplex NArray with
    `shape = [.., <%if d>2%>nz, <%end;if d>1%>ny, <%end%>nx/2+1]`
  @raise  [Numo::FFTE::RadixError] if `nx`, `ny`<% if d>2%>, `nz`<%end%>
    is not `(2^p)*(3^q)*(5^r)`.
*/
<% $funcs.push func="dzfft#{d}d" %>
static VALUE
numo_ffte_<%=func%>(int argc, VALUE *args, VALUE mod)
{
    narray_t *na;
    VALUE viopt;
    VALUE vb, vna;
    VALUE vres;
    int ndim;
    integer iopt=-1, iopt_zero=0;
    fft_opt_t *g;
    integer <%=argmap(d){|i|"n#{i}"}%>;
    size_t n=1;
    size_t shape[<%=d%>];
    ndfunc_arg_in_t ain[1] = {{cDF,<%=d%>}};
    ndfunc_arg_out_t aout[1] = {{cDC,<%=d%>,shape}};
    ndfunc_t ndf = { iter_fft_dzfft<%=d%>d, NO_LOOP, 1, 1, ain, aout };

    switch(rb_scan_args(argc, args, "11", &vna, &viopt)) {
    case 2:
        iopt = NUM2INT(viopt);
    }
    GetNArray(vna,na);
    ndim = NA_NDIM(na);
    if (ndim<<%=d%>) {
        rb_raise(eRadixError,"ndim(=%d) should >= <%=d%>",ndim);
    }
<% (1..d).each do |i| %>
    n<%=i%> = NA_SHAPE(na)[NA_NDIM(na)-<%=i%>];
    n *= n<%=i%>;
<% if i==1 %>
    shape[<%=d-1%>] = n<%=i%>/2+1;
<% else %>
    shape[<%=d-i%>] = n<%=i%>;
<% end %>
    if (!is235radix(n<%=i%>)) {
        rb_raise(eRadixError,"%d-th dim length(=%ld) is not 2,3,5-radix",
                 NA_NDIM(na)-<%=i%>,shape[<%=d-i%>]);
    }
<% end %>

    vna = na_copy(vna);
    g = alloc_fft_opt(n, iopt, &vb);

    dzfft<%=d%>d_(NULL, <%=argmap(d){|i|"&n#{i}"}%>, &iopt_zero, g->buf);
    vres = na_ndloop3(&ndf, g, 1, vna);
    RB_GC_GUARD(vb);

    return vres;
}
<% end %>


void
Init_ffte()
{
    VALUE mNumo,mFFTE;

    rb_require("numo/narray");
    mNumo = rb_define_module("Numo");
    mFFTE = rb_define_module_under(mNumo,"FFTE");
    // Radix Error
    eRadixError = rb_define_class_under(mFFTE,"RadixError",rb_eStandardError);

<% $funcs.each do |f| %>
    rb_define_module_function(mFFTE,"<%=f%>",numo_ffte_<%=f%>,-1);
<% end %>
}
