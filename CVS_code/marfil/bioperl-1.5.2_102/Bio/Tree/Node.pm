# $Id: Node.pm,v 1.49.4.3 2006/11/27 15:46:48 sendu Exp $
#
# BioPerl module for Bio::Tree::Node
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Node - A Simple Tree Node

=head1 SYNOPSIS

    use Bio::Tree::Node;
    my $nodeA = new Bio::Tree::Node();
    my $nodeL = new Bio::Tree::Node();
    my $nodeR = new Bio::Tree::Node();

    my $node = new Bio::Tree::Node();
    $node->add_Descendent($nodeL);
    $node->add_Descendent($nodeR);

    print "node is not a leaf \n" if( $node->is_leaf);

=head1 DESCRIPTION

Makes a Tree Node suitable for building a Tree.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Aaron Mackey, amackey-at-virginia-dot-edu
Sendu Bala,   bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::Node;
use vars qw($CREATIONORDER);
use strict;

use Scalar::Util qw(weaken isweak);

use base qw(Bio::Root::Root Bio::Tree::NodeI);

BEGIN { 
    $CREATIONORDER = 0;
}

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::Node();
 Function: Builds a new Bio::Tree::Node object
 Returns : Bio::Tree::Node
 Args    : -descendents   => arrayref of descendents (they will be
                             updated s.t. their ancestor point is this
                             node)
           -branch_length => branch length [integer] (optional)
           -bootstrap     => value   bootstrap value (string)
           -description   => description of node
           -id            => human readable id for node

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($children, $branchlen,$id,
      $bootstrap, $desc,$d) = $self->_rearrange([qw(
						    DESCENDENTS
						    BRANCH_LENGTH
						    ID
						    BOOTSTRAP
						    DESC
						    DESCRIPTION
						 )],
					     @args);
  $self->_register_for_cleanup(\&node_cleanup);
  $self->{'_desc'} = {}; # for descendents
  if( defined $d && defined $desc ) { 
      $self->warn("can only accept -desc or -description, not both, accepting -description");
      $desc = $d;
  } elsif( defined $d && ! defined $desc ) {
      $desc = $d;
  }
  defined $desc && $self->description($desc);
  defined $bootstrap && $self->bootstrap($bootstrap);
  defined $id && $self->id($id);
  defined $branchlen && $self->branch_length($branchlen);
  if( defined $children ) {
      if( ref($children) !~ /ARRAY/i ) {
	  $self->warn("Must specify a valid ARRAY reference to initialize a Node's Descendents");
      }
      foreach my $c ( @$children ) { 	
	  $self->add_Descendent($c);
      }
  }
  $self->_creation_id($CREATIONORDER++);
  return $self;
}

=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $node->add_Descendant($node);
 Function: Adds a descendent to a node
 Returns : number of current descendents for this node
 Args    : Bio::Node::NodeI
           boolean flag, true if you want to ignore the fact that you are
           adding a second node with the same unique id (typically memory 
           location reference in this implementation).  default is false and 
           will throw an error if you try and overwrite an existing node.

=cut

sub add_Descendent{
   my ($self,$node,$ignoreoverwrite) = @_;
   return -1 if( ! defined $node );
   
   if( ! ref($node) ||
       ref($node) =~ /HASH/ ||
       ! $node->isa('Bio::Tree::NodeI') ) {
       $self->throw("Trying to add a Descendent who is not a Bio::Tree::NodeI");
       return -1;
   }
   
   $self->{_adding_descendent} = 1;
   # avoid infinite recurse
   $node->ancestor($self) unless $node->{_setting_ancestor}; 
   $self->{_adding_descendent} = 0;
   
   if( $self->{'_desc'}->{$node->internal_id} && ! $ignoreoverwrite ) {
       $self->throw("Going to overwrite a node which is $node that is already stored here, set the ignore overwrite flag (parameter 2) to true to ignore this in the future");
   }
   $self->{'_desc'}->{$node->internal_id} = $node; # is this safely unique - we've tested before at any rate??
   
   $self->invalidate_height();
   
   return scalar keys %{$self->{'_desc'}};
}

=head2 each_Descendent

 Title   : each_Descendent($sortby)
 Usage   : my @nodes = $node->each_Descendent;
 Function: all the descendents for this Node (but not their descendents
					      i.e. not a recursive fetchall)
 Returns : Array of Bio::Tree::NodeI objects
 Args    : $sortby [optional] "height", "creation", "alpha", "revalpha",
           or coderef to be used to sort the order of children nodes.

=cut

sub each_Descendent{
   my ($self, $sortby) = @_;

   # order can be based on branch length (and sub branchlength)   
   $sortby ||= 'none';
   if (ref $sortby eq 'CODE') {
       my @values = sort { $sortby->($a,$b) } values %{$self->{'_desc'}};
       return @values;
   } elsif ($sortby eq 'height') {
       return map { $_->[0] }
       sort { $a->[1] <=> $b->[1] || 
		  $a->[2] <=> $b->[2] } 
       map { [$_, $_->height, $_->internal_id ] } 
       values %{$self->{'_desc'}};
   } elsif( $sortby eq 'alpha' ) {
       my @set;
       for my $v ( values %{$self->{'_desc'}} ) {
	   unless( $v->is_Leaf ) {
	       my @lst = ( sort { $a cmp $b } map { $_->id } 
                          grep { $_->is_Leaf } 
			   $v->get_all_Descendents($sortby));
	       push @set, [$v, $lst[0], $v->internal_id];
	   } else {
	       push @set, [$v, $v->id, $v->internal_id];
	   }
       } 
       return map { $_->[0] }
       sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @set;       
   } elsif( $sortby eq 'revalpha' ) {
       my @set;
       for my $v ( values %{$self->{'_desc'}} ) {
	   if( ! defined $v->id && 
	       ! $v->is_Leaf ) {
	       my ($l) = ( sort { $b cmp $a } map { $_->id } 
			   grep { $_->is_Leaf } 
			   $v->get_all_Descendents($sortby));
	       push @set, [$v, $l, $v->internal_id];
	   } else { 
	       push @set, [$v, $v->id, $v->internal_id];
	   }
       } 
       return map { $_->[0] }
       sort {$b->[1] cmp $a->[1] || $b->[2] <=> $a->[2] } @set;
   } else { # creation
       return map { $_->[0] }
       sort { $a->[1] <=> $b->[1] } 
       map { [$_, $_->internal_id ] }
       values %{$self->{'_desc'}};	   
   }
}

=head2 remove_Descendent

 Title   : remove_Descendent
 Usage   : $node->remove_Descedent($node_foo);
 Function: Removes a specific node from being a Descendent of this node
 Returns : nothing
 Args    : An array of Bio::Node::NodeI objects which have been previously
           passed to the add_Descendent call of this object.

=cut

sub remove_Descendent{
   my ($self,@nodes) = @_;
   my $c= 0;
   foreach my $n ( @nodes ) { 
       if( $self->{'_desc'}->{$n->internal_id} ) {
        $self->{_removing_descendent} = 1;
        $n->ancestor(undef);
        $self->{_removing_descendent} = 0;
	   # should be redundant
	   $self->{'_desc'}->{$n->internal_id}->ancestor(undef);
	   delete $self->{'_desc'}->{$n->internal_id};
	   $c++;
       } else { 
	   if( $self->verbose ) {
	       $self->debug(sprintf("no node %s (%s) listed as a descendent in this node %s (%s)\n",$n->id, $n,$self->id,$self));
	       $self->debug("Descendents are " . join(',', keys %{$self->{'_desc'}})."\n");
	   }
       }
   }
   
   # remove unecessary nodes if we have removed the part 
   # which branches.
   my $a1 = $self->ancestor;   
   if( $a1 ) {
       my $bl = $self->branch_length || 0;
       my @d = $self->each_Descendent;
       if (scalar @d == 1) {
	   $d[0]->branch_length($bl + ($d[0]->branch_length || 0));
	   $a1->add_Descendent($d[0]);
	   $a1->remove_Descendent($self);
       }
   }
   $c;
}

=head2 remove_all_Descendents

 Title   : remove_all_Descendents
 Usage   : $node->remove_All_Descendents()
 Function: Cleanup the node's reference to descendents and reset
           their ancestor pointers to undef, if you don't have a reference
           to these objects after this call they will be cleaned up - so
           a get_nodes from the Tree object would be a safe thing to do first
 Returns : nothing
 Args    : none

=cut

sub remove_all_Descendents{
   my ($self) = @_;
   # this won't cleanup the nodes themselves if you also have
   # a copy/pointer of them (I think)...
   while( my ($node,$val) = each %{ $self->{'_desc'} } ) {
       $val->ancestor(undef);
   }
   $self->{'_desc'} = {};
   1;
}

=head2 get_all_Descendents

 Title   : get_all_Descendents
 Usage   : my @nodes = $node->get_all_Descendents;
 Function: Recursively fetch all the nodes and their descendents
           *NOTE* This is different from each_Descendent
 Returns : Array or Bio::Tree::NodeI objects
 Args    : none

=cut

# get_all_Descendents implemented in the interface 

=head2 ancestor

 Title   : ancestor
 Usage   : $obj->ancestor($newval)
 Function: Set the Ancestor
 Returns : ancestral node
 Args    : newvalue (optional)

=cut

sub ancestor {
    my $self = shift;
    if (@_) {
        my $new_ancestor = shift;
        
        # we can set ancestor to undef
        if ($new_ancestor) {
            $self->throw("This is [$new_ancestor], not a Bio::Tree::NodeI")
		unless $new_ancestor->isa('Bio::Tree::NodeI');
        }
        
        my $old_ancestor = $self->{'_ancestor'} || '';
        if (!$old_ancestor || 
	    ($old_ancestor && ( !$new_ancestor || 
			       $new_ancestor ne $old_ancestor)) ) {
            if( $old_ancestor && ! $old_ancestor->{_removing_descendent}) {
		$old_ancestor->remove_Descendent($self);
	    }
            if ($new_ancestor && 
		! $new_ancestor->{_adding_descendent} ) { # avoid infinite recurse
                $self->{_setting_ancestor} = 1;
                $new_ancestor->add_Descendent($self, 1);
                $self->{_setting_ancestor} = 0;
            }
        }
        weaken($self->{'_ancestor'} = $new_ancestor);
    }
    
    return $self->{'_ancestor'};
}

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length()
 Function: Get/Set the branch length
 Returns : value of branch_length
 Args    : newvalue (optional)

=cut

sub branch_length{
    my $self = shift;
    if( @_ ) {
	my $bl = shift;
	if( defined $bl &&
	    $bl =~ s/\[(\d+)\]// ) {
	    $self->bootstrap($1);
	}
	$self->{'_branch_length'} = $bl;
    $self->invalidate_height();
    }
    return $self->{'_branch_length'};
}

=head2 bootstrap

 Title   : bootstrap
 Usage   : $obj->bootstrap($newval)
 Function: Get/Set the bootstrap value
 Returns : value of bootstrap
 Args    : newvalue (optional)

=cut

sub bootstrap { 
    my $self = shift;
    if( @_ ) {
	if( $self->has_tag('B') ) {
	    $self->remove_tag('B');
	}
	$self->add_tag_value('B',shift);
    }
    return ($self->get_tag_values('B'))[0];
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: Get/Set the description string
 Returns : value of description
 Args    : newvalue (optional)

=cut

sub description {
    my $self = shift;
    $self->{'_description'} = shift @_ if @_;
    return $self->{'_description'};
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: The human readable identifier for the node 
 Returns : value of human readable id
 Args    : newvalue (optional)

"A name can be any string of printable characters except blanks,
colons, semicolons, parentheses, and square brackets. Because you may
want to include a blank in a name, it is assumed that an underscore
character ("_") stands for a blank; any of these in a name will be
converted to a blank when it is read in."  

from L<http://evolution.genetics.washington.edu/phylip/newicktree.html>

Also note that these objects now support spaces, ();: because we can
automatically quote the strings if they contain these characters.  The
L<id_output> method does this for you so use the id() method to get
the raw string while L<id_output> to get the pre-escaped string.

=cut

sub id {
    my ($self, $value) = @_;
    if (defined $value) {
        #$self->warn("Illegal characters ();:  and space in the id [$value], converting to _ ")
	# if $value =~ /\(\);:/ and $self->verbose >= 0;
        #$value =~ s/[\(\);:\s]/_/g;
        $self->{'_id'} = $value;
    }
    return $self->{'_id'};
}

=head2 Helper Functions

=cut

=head2 id_output

 Title   : id_output
 Usage   : my $id = $node->id_output;
 Function: Return an id suitable for output in format like newick
           so that if it contains spaces or ():; characters it is properly 
           quoted
 Returns : $id string if $node->id has a value
 Args    : none

=cut

# implemented in NodeI interface 

=head2 internal_id

 Title   : internal_id
 Usage   : my $internalid = $node->internal_id
 Function: Returns the internal unique id for this Node
           (a monotonically increasing number for this in-memory implementation
            but could be a database determined unique id in other 
	    implementations)
 Returns : unique id
 Args    : none

=cut

sub internal_id {
   return $_[0]->_creation_id;
}

=head2 _creation_id

 Title   : _creation_id
 Usage   : $obj->_creation_id($newval)
 Function: a private method signifying the internal creation order
 Returns : value of _creation_id
 Args    : newvalue (optional)

=cut

sub _creation_id {
    my $self = shift @_;
    $self->{'_creation_id'} = shift @_ if( @_);
    return $self->{'_creation_id'} || 0;
}

=head2 Bio::Node::NodeI decorated interface implemented

The following methods are implemented by L<Bio::Node::NodeI> decorated
interface.

=head2 is_Leaf

 Title   : is_Leaf
 Usage   : if( $node->is_Leaf )
 Function: Get Leaf status
 Returns : boolean
 Args    : none

=cut

sub is_Leaf {
    my ($self) = @_;
    my $isleaf = ! (defined $self->{'_desc'} &&
		 (keys %{$self->{'_desc'}} > 0) );
    return $isleaf;
}

=head2 to_string

 Title   : to_string
 Usage   : my $str = $node->to_string()
 Function: For debugging, provide a node as a string
 Returns : string
 Args    : none

=head2 height

 Title   : height
 Usage   : my $len = $node->height
 Function: Returns the height of the tree starting at this
           node.  Height is the maximum branchlength to get to the tip.
 Returns : The longest length (weighting branches with branch_length) to a leaf
 Args    : none

=cut

sub height { 
    my ($self) = @_;
    return $self->{'_height'} if( defined $self->{'_height'} );
    
    return 0 if( $self->is_Leaf );
    my $max = 0;
    foreach my $subnode ( $self->each_Descendent ) { 
	my $bl = $subnode->branch_length;
	$bl = 1 unless (defined $bl && $bl =~ /^\-?\d+(\.\d+)?$/);
	my $s = $subnode->height + $bl;
	if( $s > $max ) { $max = $s; }
    }
    return ($self->{'_height'} = $max);
}

=head2 invalidate_height

 Title   : invalidate_height
 Usage   : private helper method
 Function: Invalidate our cached value of the node height in the tree
 Returns : nothing
 Args    : none

=cut

sub invalidate_height { 
    my ($self) = @_;
    
    $self->{'_height'} = undef;
    if( defined $self->ancestor ) {
	$self->ancestor->invalidate_height;
    }
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $node->add_tag_value($tag,$value)
 Function: Adds a tag value to a node 
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub add_tag_value{
    my ($self,$tag,$value) = @_;
    if( ! defined $tag || ! defined $value ) {
	$self->warn("cannot call add_tag_value with an undefined value");
    }
    push @{$self->{'_tags'}->{$tag}}, $value;
    return scalar @{$self->{'_tags'}->{$tag}};
}

=head2 remove_tag

 Title   : remove_tag
 Usage   : $node->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove


=cut

sub remove_tag {
   my ($self,$tag) = @_;
   if( exists $self->{'_tags'}->{$tag} ) {
       $self->{'_tags'}->{$tag} = undef;
       delete $self->{'_tags'}->{$tag};
       return 1;
   }
   return 0;
}

=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $node->remove_all_tags()
 Function: Removes all tags 
 Returns : None
 Args    : None

=cut

sub remove_all_tags{
   my ($self) = @_;
   $self->{'_tags'} = {};
   return;
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $node->get_all_tags()
 Function: Gets all the tag names for this Node
 Returns : Array of tagnames
 Args    : None

=cut

sub get_all_tags{
   my ($self) = @_;
   my @tags = sort keys %{$self->{'_tags'} || {}};
   return @tags;
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $node->get_tag_value($tag)
 Function: Gets the values for given tag ($tag)
 Returns : Array of values or empty list if tag does not exist
 Args    : $tag - tag name

=cut

sub get_tag_values{
   my ($self,$tag) = @_;
   return wantarray ? @{$self->{'_tags'}->{$tag} || []} :
                     (@{$self->{'_tags'}->{$tag} || []})[0];
}

=head2 has_tag

 Title   : has_tag
 Usage   : $node->has_tag($tag)
 Function: Boolean test if tag exists in the Node
 Returns : Boolean
 Args    : $tag - tagname

=cut

sub has_tag {
   my ($self,$tag) = @_;
   return exists $self->{'_tags'}->{$tag};
}

sub node_cleanup {
    my $self = shift;
    return unless defined $self;
    
    #*** below is wrong, cleanup doesn't actually occur. Will replace with:
    # $self->remove_all_Descendents; once further fixes in place..
    if( defined $self->{'_desc'} &&
        ref($self->{'_desc'}) =~ /ARRAY/i ) {
        while( my ($nodeid,$node) = each %{ $self->{'_desc'} } ) {
            $node->ancestor(undef); # insure no circular references
            $node = undef;
        }
    }
    $self->{'_desc'} = {};
}

=head2 reverse_edge

 Title   : reverse_edge
 Usage   : $node->reverse_edge(child);
 Function: makes child be a parent of node
 Requires: child must be a direct descendent of node
 Returns : 1 on success, 0 on failure
 Args    : Bio::Tree::NodeI that is in the tree

=cut

sub reverse_edge {
    my ($self,$node) = @_;
    if( $self->delete_edge($node) ) {
      $node->add_Descendent($self);
      return 1;
    } 
    return 0;
}

1;
