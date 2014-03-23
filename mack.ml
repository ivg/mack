open Batteries
module M = Gsl.Matrix
module V = Gsl.Vector
module A = Array
module BLAS = Gsl.Blas
module Big = Bigarray
module BA = Big.Array2

open Format

type problem = Maxi | Mini
type matrix = M.matrix

let normalized_copy _M problem =
  let worst,scaler = match problem with
    | Maxi -> neg_infinity, ~-.1.
    | Mini -> infinity, 1. in
  let xform = function
    | e when e <> worst -> e *. scaler
    | _ -> max_float in
  BA.map xform Big.float64 _M


let corrected_copy _M problem =
  let _M = M.copy _M in
  match problem with
    | Maxi -> M.scale _M ~-.1.; _M
    | Mini -> _M



module ISet :
sig
  type t
  val empty: int -> t
  val singleton: int -> int -> t
  val mem: int -> t -> bool
  val is_empty: t -> bool
  val is_singleton: t -> bool
  val size: t -> int
  val add: int -> t -> t
  val rem: int -> t -> t
  val iter: (int -> unit) -> t -> unit
  val fold_left: ('a -> int -> 'a) -> 'a -> t -> 'a
  val union: t -> t -> t
  val compl: t -> t
  val first: t -> int
  val fold: ('a -> int -> 'a) -> 'a -> t -> 'a
end =
struct
  type ilist = int list
  type bitset = BitSet.t
  type t = ilist * ilist * bitset * bitset
  let empty n =
    [],List.init n identity,BitSet.create n,BitSet.create_full n
  let mem i (_,_,s,_) = BitSet.mem s i
  let is_empty = function ([],_,_,_) -> true | _ -> false
  let is_singleton = function (_::[],_,_,_) -> true | _ -> false
  let size (n,_,_,_) = List.length n
  let add i (n,n',s,s') =
    let s,s' = BitSet.copy s, BitSet.copy s' in
    BitSet.set s i;
    BitSet.unset s' i;
    i::n,List.remove n' i, s, s'

  let rem i (n,n',s,s') =
    let s,s' = BitSet.copy s, BitSet.copy s' in
    BitSet.unset s i;
    BitSet.set s' i;
    List.remove n i, i::n', s, s'

  let iter f (n,_,_,_) = List.iter f n
  let fold_left f z (n,_,_,_) = List.fold_left f z n
  let singleton n e = let s = empty n in add e s

  let union (ln,ln',ls,ls') (rn,rn',rs,rs') =
    (* XXX: объединение намеренно сломано в целях оптимизации
       оно работает правильно только для Мака.
    *)
    let s  = BitSet.union ls rs in
    let s' = BitSet.union ls' rs' in
    let n  = ln @ rn  in
    let n' = [-1] in
    n,n',s,s'

  let compl (n,n',s,s') = (n',n,s',s)
  let first (n,_,_,_) = List.hd n
  let fold f init (n,_,_,_) = List.fold_left f init n
end

module RSet = struct
  type t = ISet.t array
  let create n = A.init n (fun _ -> ISet.empty n)
  let add i j s = s.(j) <- ISet.add i s.(j)
  let rem i j s = s.(j) <- ISet.rem i s.(j)
  let rows j s = s.(j)
  let all_singleton s = A.for_all ISet.is_singleton s
  let find_nonsingleton s =
    A.findi (fun r -> not (ISet.is_singleton r || ISet.is_empty r)) s
end


let solve prob _M =
  let _M = normalized_copy _M prob in
  let m,_ = M.dims _M in
  let last = m - 1 in
  let base = A.make m ~-1 in
  let rows = RSet.create m in

  let select_bases () =
    for i = 0 to last do
      let j = V.min_index (M.row _M i) in
      base.(i) <- j;
      RSet.add i j rows
    done in

  let have_solution () = RSet.all_singleton rows in

  let find_min _ROW _COL =
    let _COL = ISet.compl _COL in
    ISet.fold (fun s i -> ISet.fold (fun (rr,kk,d) j ->
      let r = _M.{i,j} -. _M.{i, base.(i)} in
      if r < d then (i,j,r) else (rr,kk,d))
      s _COL) (-1, -1, infinity) _ROW in

  let adjust_dual_solution _COL delta =
    ISet.iter (fun j ->
      for i = 0 to last do
        _M.{i,j} <- _M.{i,j} +. delta
      done
    ) _COL in

  let switch_columns rr nk =
    let ok = base.(rr) in
    base.(rr) <- nk;
    RSet.rem rr ok rows;
    RSet.add rr nk rows in

  let rec main_loop () =
    if have_solution () then base
    else
      let j = RSet.find_nonsingleton rows in
      let _COL = ISet.singleton m j and _ROW = RSet.rows j rows in
      inner_loop _ROW _COL
  and inner_loop _ROW _COL =
    let rr,kk,delta = find_min _ROW _COL in
    adjust_dual_solution _COL delta;
    match ISet.is_empty (RSet.rows kk rows) with
      | true -> begin
        switch_columns rr kk;
        main_loop ()
      end
      | false  -> begin
        switch_columns rr kk;
        let _COL = ISet.add kk _COL in
        let _ROW = ISet.union (RSet.rows kk rows) _ROW in
        inner_loop _ROW _COL
      end in
  select_bases ();
  main_loop ()



