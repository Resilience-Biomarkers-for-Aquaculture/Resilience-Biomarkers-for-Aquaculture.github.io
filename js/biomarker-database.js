$(document).ready(function() {
  var filterColIndices = [2, 3, 4, 5, 6, 7]; // species, condition, trait, effect, class, evidence

  var table = $('#biomarker-table').DataTable({
    dom: '<"row"<"col-sm-4"B><"col-sm-4"l><"col-sm-4"f>>' +
         '<"row"<"col-sm-12"tr>>' +
         '<"row"<"col-sm-5"i><"col-sm-7"p>>',
    buttons: [
      {
        extend: 'csvHtml5',
        text: '<span class="glyphicon glyphicon-download-alt"></span> CSV',
        title: 'biomarker_database',
        className: 'btn-sm btn-default'
      },
      {
        extend: 'excelHtml5',
        text: '<span class="glyphicon glyphicon-download-alt"></span> Excel',
        title: 'biomarker_database',
        className: 'btn-sm btn-default'
      },
      {
        extend: 'copyHtml5',
        text: '<span class="glyphicon glyphicon-copy"></span> Copy',
        className: 'btn-sm btn-default'
      }
    ],
    pageLength: 25,
    orderCellsTop: true,
    order: [[9, 'desc']],
    columnDefs: [
      { targets: 9, type: 'num' },
      {
        // Linkify URLs in the Source column for display; keep plain text for search/sort.
        targets: 8,
        render: function(data, type) {
          if (type !== 'display') { return data; }
          return data.replace(
            /\((https?:\/\/[^\s)]+)\)/g,
            function(match, url) {
              // Guard: only allow http/https URLs to prevent javascript: or other protocol injection.
              if (!/^https?:\/\//i.test(url)) { return match; }
              var safeHref = encodeURI(url);
              var safeText = url.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
              return '(<a href="' + safeHref + '" target="_blank" rel="noopener noreferrer">' + safeText + '</a>)';
            }
          );
        }
      }
    ],
    initComplete: function() {
      var api = this.api();
      var filterRow = $('#biomarker-table thead tr:eq(1)');
      // Per-column dropdown filters for categorical columns.
      filterColIndices.forEach(function(colIdx) {
        var column = api.column(colIdx);
        var headerCell = filterRow.find('th:eq(' + colIdx + ')');
        var select = $('<select class="form-control input-sm"><option value="">-- All --</option></select>')
          .appendTo($(headerCell).empty())
          .on('change', function() {
            var val = $.fn.dataTable.util.escapeRegex($(this).val());
            column.search(val ? '^' + val + '$' : '', true, false).draw();
          });
        column.data().unique().sort().each(function(d) {
          if (d) { select.append($('<option>').val(d).text(d)); }
        });
      });
    }
  });
});
